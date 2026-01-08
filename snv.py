"""
TreeHarmonizer SNV Placement Module

This module handles placement of single nucleotide variants (SNVs) onto the
phylogenetic tree, with optional regenotyping using structural variants (SVs)
and/or copy number alterations (CNAs).

Workflow Overview:
    1. Load SNV data and phylogenetic tree
    2. Calculate clade-based support requirements using false negative rate
    3. (Optional) Process SV/CNA data for regenotyping evidence
    4. Apply regenotyping logic to expand variant sublines
    5. Perform category scoring to determine final placements
    6. Generate exclusive and cumulative VCF files per tree node

Regenotyping Modes:
    - 'severus': Use only SV deletion data for regenotyping
    - 'cna': Use only CNA (CN0/CN1) data for regenotyping
    - 'both': Use both SV and CNA data (default)

All modes implement "dramatic shift" prevention to avoid placing variants
at nodes with >2x the original subline count.

Key Concepts:
    - DV_SUBLINES: Samples where the variant was called by DeepVariant
    - MRCA: Most Recent Common Ancestor - the tree node where variant is placed
    - Parsimony Assumption: SV/CNA must NOT co-occur with SNV on same haplotype
      (if they do, regenotyping doesn't apply because deletion removes the SNV)
    - Support Threshold: Minimum samples needed based on clade size and FN rate
"""

import pandas as pd
import subprocess
import utils as th_utils


# =============================================================================
# MAIN ENTRY POINT
# =============================================================================

def snv_placement_runner(snv_path, tree_newick, tree_metadata=None, output_path=None, fn_rate=0.15,
                         do_regenotyping=False, sv_path=None, cna_path=None,
                         regenotype_mode='both',
                         write_exclusive_vcfs=True, write_cumulative_vcfs=True,
                         sample_list=None, chromosomes=None, verbose=False,
                         command_string=None, version=None):
    """
    Run SNV placement with optional regenotyping using SV and/or CNA data.

    This is the main entry point for SNV placement. It orchestrates loading data,
    calculating support requirements, optionally applying regenotyping, and
    generating output VCF files.

    Args:
        snv_path (str): Path to directory containing per-sample DeepVariant VCF files.
        tree_newick (str): Newick-format string representing the phylogenetic tree.
        tree_metadata (TreeMetadata): Metadata about tree structure and node mappings.
            If None, will be auto-generated from tree_newick using normalize_tree().
        output_path (str): Base output directory for results and VCFs.
        fn_rate (float): False negative rate for support threshold calculation.
            Default 0.15 means we expect 15% of true variants to be missed,
            so we require (1 - fn_rate) * clade_size samples for placement.
        do_regenotyping (bool): If True, use SV/CNA data to expand variant sublines.
        sv_path (str): Path to Severus VCF file. Required if regenotyping with SVs.
        cna_path (str): Path to Wakhan CNA output directory. Required if regenotyping
            with CNAs. Should contain per-sample subdirectories with BED files.
        regenotype_mode (str): Which data sources to use for regenotyping:
            - 'severus': SV deletions only
            - 'cna': CNA (CN0/CN1 regions) only
            - 'both': Both SV and CNA data
        write_exclusive_vcfs (bool): If True, write VCFs with variants exclusively
            placed at each tree node.
        write_cumulative_vcfs (bool): If True, write VCFs with variants at each
            node plus all descendant nodes.
        sample_list (list): Optional list of sample names. If None, auto-detected
            from input files.
        chromosomes (list): Optional list of chromosomes to process.
        verbose (bool): If True, print detailed per-node variant counts.
        command_string (str): Optional command line string for provenance tracking
            in output VCF headers.

    Returns:
        dict: Results dictionary containing:
            - 'marked_df': DataFrame with all variants and placement metadata
            - 'exclusive_dfs': Dict mapping node names to exclusive variant DataFrames
            - 'cumulative_dfs': Dict mapping node names to cumulative variant DataFrames
            - 'placed_count': Number of variants successfully placed
            - 'unplaced_count': Number of variants that failed support threshold

    Notes:
        - Writes placed/unplaced variant TSV files to {output_path}/snv/
        - Writes exclusive VCFs to {output_path}/snv/exclusive/
        - Writes cumulative VCFs to {output_path}/snv/cumulative/
    """
    from tree_preprocessing import normalize_tree

    # Auto-generate tree_metadata if not provided (backward compatibility)
    if tree_metadata is None:
        tree_newick, tree_metadata = normalize_tree(tree_newick, verbose=verbose)

    print("Loading SNV data and tree...")

    # Configure pandas display options
    pd.set_option('display.width', 5000)
    pd.set_option("display.expand_frame_repr", True)
    pd.set_option("display.max_colwidth", 1000)
    pd.set_option('display.max_columns', None)

    # Load SNV data into a merged DataFrame
    dv_merged, detected_sample_list, vcf_metadata = th_utils.generate_merged_df(
        caller_path=snv_path,
        predefined_sample_list=sample_list,
        chromosomes=chromosomes
    )

    # Print any FORMAT warnings if verbose
    if verbose and vcf_metadata.get('format_warnings'):
        print("Warning: FORMAT field inconsistencies detected:")
        for warning in vcf_metadata['format_warnings']:
            print(f"  - {warning}")

    # Use detected samples if not provided
    if sample_list is None:
        sample_list = detected_sample_list

    # Load tree input via newick string and parse it into various components
    imported_tree, _, _, non_terminal_paths, _, _, _, _ = \
        th_utils.get_tree_data(tree_newick, tree_metadata)

    # Add informative columns that were lost from original merging
    dv_merged['CHROM'] = dv_merged['KEY'].str.split(":").str[0]
    dv_merged['POS'] = dv_merged['KEY'].str.split(":").str[1]
    dv_merged['REF'] = dv_merged['KEY'].str.split(":").str[2]
    dv_merged['ALT'] = dv_merged['KEY'].str.split(":").str[3]

    dv_merged['CHROM'] = dv_merged['CHROM'].astype(str)
    dv_merged['POS'] = dv_merged['POS'].astype(int)

    # Rename INFO column to avoid conflicts
    dv_merged.rename(columns={'INFO': 'DV_INFO'}, inplace=True)

    print(f"Loaded {len(dv_merged)} SNV variants across {len(sample_list)} samples")

    # Initialize MARKED_dv which will contain all metadata
    MARKED_dv = dv_merged.copy(deep=True)

    # Calculate DV_SUBLINES early - needed for both regenotyping and non-regenotyping paths
    MARKED_dv['DV_SUBLINES'] = MARKED_dv.apply(
        lambda row: [col for col in sample_list if pd.notna(row[col])],
        axis=1
    )

    # Calculate clade size support requirements
    # This determines how many samples must have a variant for it to be placed at each clade.
    # The formula accounts for expected false negatives in variant calling.
    print("Calculating clade size support requirements...")
    clade_sizes = set()
    for clade in imported_tree.traverse():
        if clade.is_leaf():
            continue
        clade_sizes.add(len(clade.get_leaves()))

    # Build lookup table: clade_size -> minimum samples required
    # Formula: For clade size N, require ceil(N * (1 - fn_rate)) samples
    # Special cases: size 1 requires 1, size 2 requires 2 (no tolerance for FN)
    minimum_subline_support_per_clade_size_requirement = {}
    for clade_size in clade_sizes:
        if clade_size < 2:
            minimum_subline_support_per_clade_size_requirement[clade_size] = 1
        elif clade_size == 2:
            # Size-2 clades require both samples (stringent to avoid false placements)
            minimum_subline_support_per_clade_size_requirement[clade_size] = 2
        else:
            # For larger clades, allow some missing due to false negatives
            # e.g., with fn_rate=0.15 and clade_size=10, require 8 samples
            support_requirement = int(clade_size * (1 - float(fn_rate)))
            minimum_subline_support_per_clade_size_requirement[clade_size] = support_requirement

    # Process regenotyping if enabled
    if do_regenotyping:
        print(f"\n=== Regenotyping Mode: {regenotype_mode.upper()} ===")

        # Determine which data sources are needed
        need_sv = regenotype_mode in ['severus', 'both']
        need_cna = regenotype_mode in ['cna', 'both']

        # Process SV data if needed
        if need_sv and sv_path:
            print("Processing Severus SV data for regenotyping...")
            MARKED_dv = process_severus_for_regenotyping(
                MARKED_dv, sv_path, sample_list, imported_tree, tree_metadata, chromosomes=chromosomes
            )

        # Process CNA data if needed
        if need_cna and cna_path:
            print("Processing Wakhan CNA data for regenotyping...")
            MARKED_dv = process_wakhan_for_regenotyping(
                MARKED_dv, cna_path, sample_list, imported_tree, tree_metadata, chromosomes=chromosomes
            )

        # Apply regenotyping logic
        print("Applying regenotyping logic...")
        MARKED_dv = apply_regenotyping_logic(
            MARKED_dv, regenotype_mode, minimum_subline_support_per_clade_size_requirement,
            sample_list, imported_tree, tree_metadata
        )
    else:
        # Just calculate original DV placement (DV_SUBLINES already calculated above)
        print("Calculating original SNV placement (no regenotyping)...")
        MARKED_dv['DV_MRCA'] = MARKED_dv.apply(
            lambda row: th_utils.common_ancestor_helper(row, "DV_SUBLINES", input_tree=imported_tree, tree_metadata=tree_metadata),
            axis=1
        )
        MARKED_dv['DV_MRCA_TERMINALS'] = MARKED_dv.apply(
            lambda row: imported_tree.search_nodes(name=row['DV_MRCA'])[0].get_leaf_names(),
            axis=1
        )
        MARKED_dv['DV_MINIMUM_SUPPORT_MET'] = MARKED_dv.apply(
            lambda row: len(row['DV_SUBLINES']) >=
                minimum_subline_support_per_clade_size_requirement[len(row['DV_MRCA_TERMINALS'])],
            axis=1
        )
        MARKED_dv['FINAL_MRCA'] = MARKED_dv.apply(
            lambda row: row['DV_MRCA'] if row['DV_MINIMUM_SUPPORT_MET'] else float('nan'),
            axis=1
        )

    # Write debug files
    print("\nWriting analysis files...")
    subprocess.run(['mkdir', '-p', f"{output_path}/snv"])

    placed_variants = MARKED_dv[MARKED_dv['FINAL_MRCA'].notna()].copy()
    unplaced_variants = MARKED_dv[MARKED_dv['FINAL_MRCA'].isna()].copy()

    placed_variants.to_csv(f"{output_path}/snv/placed_snv_variants.tsv", sep="\t", index=False)
    unplaced_variants.to_csv(f"{output_path}/snv/unplaced_snv_variants.tsv", sep="\t", index=False)

    print(f"Placed SNV variants: {len(placed_variants)}")
    print(f"Unplaced SNV variants: {len(unplaced_variants)}")
    print(f"Placement rate: {len(placed_variants) / len(MARKED_dv) * 100:.2f}%")

    # Generate VCFs
    print("\nGenerating VCF files...")
    exclusive_dfs, cumulative_dfs = generate_vcfs(
        MARKED_dv, non_terminal_paths, output_path, sample_list,
        write_exclusive_vcfs, write_cumulative_vcfs, verbose, tree_metadata,
        vcf_metadata=vcf_metadata, command_string=command_string, version=version
    )

    #print("\nSNV placement complete!")

    return {
        'marked_df': MARKED_dv,
        'exclusive_dfs': exclusive_dfs,
        'cumulative_dfs': cumulative_dfs,
        'placed_count': len(placed_variants),
        'unplaced_count': len(unplaced_variants)
    }


# =============================================================================
# DATA PROCESSING - SEVERUS SV
# =============================================================================

def process_severus_for_regenotyping(MARKED_dv, sv_path, sample_list, imported_tree, tree_metadata, chromosomes=None):
    """
    Process Severus SV deletion data and add regenotyping metadata to MARKED_dv.

    This function loads SV deletion calls from Severus, creates interval trees
    for efficient genomic range lookups, and annotates each SNV with information
    about overlapping deletions. This metadata is used later by the regenotyping
    logic to potentially expand the sublines for each SNV.

    The key insight: If an SNV falls within a deletion region in sample X, but
    sample X doesn't have the SNV called (false negative), the deletion provides
    evidence that the SNV should be there (regenotyping).

    Args:
        MARKED_dv (pd.DataFrame): SNV DataFrame with DV_SUBLINES column already calculated.
        sv_path (str): Path to Severus merged VCF file.
        sample_list (list): List of sample names to process.
        imported_tree (ete3.Tree): Parsed phylogenetic tree object.

    Returns:
        pd.DataFrame: MARKED_dv with additional columns:
            - IN_SEVERUS_DELETION: Boolean, True if SNV overlaps any SV deletion
            - SEVERUS_SUBLINES: List of samples with the overlapping deletion
            - SEVERUS_MRCA: MRCA node of the deletion
            - SEVERUS_MRCA_TERMINALS: Leaf nodes under SEVERUS_MRCA
            - SEVERUS_ID: ID of the overlapping Severus deletion
            - SEVERUS_OTHER_HAPLO_BOOL: True if SNV and deletion share samples
              (parsimony check - see Notes)

    Notes:
        Parsimony Assumption: If a sample has BOTH the SNV called AND falls within
        the deletion region, this suggests the SNV and deletion are on different
        haplotypes. In this case, SEVERUS_OTHER_HAPLO_BOOL=True, and regenotyping
        does NOT apply for this variant because the deletion didn't remove the SNV.
    """

    # Load Severus data
    severus, sv_sample_list, _ = th_utils.generate_severus_df(
        severus_path=sv_path, simple_name=True, chromosomes=chromosomes
    )

    # Validate SV samples match SNV samples (optional but recommended)
    if set(sv_sample_list) != set(sample_list):
        print(f"Warning: SV sample list {sv_sample_list} differs from SNV sample list {sample_list}")
    severus['CHROM'] = severus['CHROM'].astype(str)
    severus['POS'] = severus['POS'].astype(int)
    severus.rename(columns={'INFO': 'SEV_INFO'}, inplace=True)

    # Filter out second breakpoint pair
    severus_no_break = severus[severus["ID"].str.contains("_2") == False]
    severus_filtered = severus_no_break.copy(deep=True)

    # Filter for only deletions
    severus_filtered_del = severus_filtered[
        severus_filtered["SEV_INFO"].str.contains("SVTYPE=DEL") == True
    ]

    # Determine called sublines for each deletion
    def severus_called_sublines_helper(row):
        output_subline_list = []
        for col in sample_list:
            internal_sev_data = row[col].split(":")
            DV = int(internal_sev_data[4])
            if DV > 0:
                output_subline_list.append(col)
        return output_subline_list

    pd.options.mode.chained_assignment = None
    severus_filtered_del['SEVERUS_SUBLINES'] = severus_filtered_del.apply(
        severus_called_sublines_helper, axis=1
    )
    severus_filtered_del['SEVERUS_MRCA'] = severus_filtered_del.apply(
        lambda row: th_utils.common_ancestor_helper(row, "SEVERUS_SUBLINES", input_tree=imported_tree, tree_metadata=tree_metadata),
        axis=1
    )
    severus_filtered_del['SEVERUS_MRCA_TERMINALS'] = severus_filtered_del.apply(
        lambda row: imported_tree.search_nodes(name=row['SEVERUS_MRCA'])[0].get_leaf_names(),
        axis=1
    )
    pd.options.mode.chained_assignment = 'warn'

    severus_filtered_del.index = severus_filtered_del['ID']

    # Create interval tree of all Severus deletion ranges for fast range queries
    # Each chromosome gets its own interval tree for efficient lookups
    # Interval trees allow O(log n + k) queries where k is number of overlapping intervals
    severus_internal_trees_per_chromosome = {}
    # Use provided chromosomes or fall back to autosomes for backward compatibility
    chrom_list = chromosomes if chromosomes else th_utils.autosomes
    for chrom in chrom_list:
        severus_internal_trees_per_chromosome[chrom] = th_utils.it.IntervalTree()

    # Populate interval trees with deletion coordinates
    # Each interval stores: (deletion_ID, MRCA_node, list_of_samples_with_deletion)
    for index, row in severus_filtered_del.iterrows():
        # Extract END position from INFO field (format: ...;END=XXXXX;...)
        end_pos = int(row['SEV_INFO'].split(";")[3].split("=")[1])
        severus_internal_trees_per_chromosome[row['CHROM']].addi(
            int(row['POS'] + 1),  # Start position (+1 for 0-based to 1-based)
            end_pos + 1,          # End position (+1 for half-open interval)
            (row['ID'], row['SEVERUS_MRCA'], row['SEVERUS_SUBLINES'])
        )

    # Add Severus metadata to MARKED_dv
    MARKED_dv['IN_SEVERUS_DELETION'] = MARKED_dv.apply(
        lambda row: len(severus_internal_trees_per_chromosome[row['CHROM']][int(row['POS'])]) > 0,
        axis=1
    )
    MARKED_dv['SEVERUS_SUBLINES'] = MARKED_dv.apply(
        lambda row: severus_internal_trees_per_chromosome[row['CHROM']][int(row['POS'])].pop()[2][2]
            if row['IN_SEVERUS_DELETION'] else float("nan"),
        axis=1
    )
    MARKED_dv['SEVERUS_MRCA'] = MARKED_dv.apply(
        lambda row: severus_internal_trees_per_chromosome[row['CHROM']][int(row['POS'])].pop()[2][1]
            if row['IN_SEVERUS_DELETION'] else float("nan"),
        axis=1
    )
    MARKED_dv['SEVERUS_MRCA_TERMINALS'] = MARKED_dv.apply(
        lambda row: imported_tree.search_nodes(name=row['SEVERUS_MRCA'])[0].get_leaf_names()
            if row['IN_SEVERUS_DELETION'] else float("nan"),
        axis=1
    )
    MARKED_dv['SEVERUS_ID'] = MARKED_dv.apply(
        lambda row: severus_internal_trees_per_chromosome[row['CHROM']][int(row['POS'])].pop()[2][0]
            if row['IN_SEVERUS_DELETION'] else float("nan"),
        axis=1
    )

    # Apply parsimony assumption check
    # If any sample has BOTH the SNV called AND the deletion, this suggests:
    # - The SNV and deletion are on DIFFERENT haplotypes (one on each chromosome copy)
    # - The deletion didn't actually remove the SNV, so regenotyping shouldn't apply
    # SEVERUS_OTHER_HAPLO_BOOL=True means regenotyping is blocked for this variant
    MARKED_dv['SEVERUS_OTHER_HAPLO_BOOL'] = MARKED_dv.apply(
        lambda row: len(set(row['DV_SUBLINES']).intersection(set(row['SEVERUS_SUBLINES']))) > 0
            if row['IN_SEVERUS_DELETION'] else float("nan"),
        axis=1
    )
    MARKED_dv['SEVERUS_OTHER_HAPLO_BOOL'] = MARKED_dv['SEVERUS_OTHER_HAPLO_BOOL'].astype('boolean')

    return MARKED_dv


# =============================================================================
# DATA PROCESSING - WAKHAN CNA
# =============================================================================

def process_wakhan_for_regenotyping(MARKED_dv, cna_path, sample_list, imported_tree, tree_metadata, chromosomes=None):
    """
    Process Wakhan CNA (copy number alteration) data and add regenotyping metadata.

    This function loads copy number segment data from Wakhan for each sample,
    creates interval trees for CN0 (homozygous deletion) and CN1 (heterozygous
    deletion/LOH) regions, and annotates SNVs with CNA overlap information.

    The key insight: If an SNV falls within a CN1 region (one copy deleted),
    a sample without the SNV call might still have the variant on the remaining
    copy - providing evidence for regenotyping.

    Args:
        MARKED_dv (pd.DataFrame): SNV DataFrame with DV_SUBLINES already calculated.
        cna_path (str): Path to Wakhan output directory containing per-sample
            subdirectories with BED files (e.g., {sample}/bed_output/{sample}_copynumbers_segments.bed).
        sample_list (list): List of sample names to process.
        imported_tree (ete3.Tree): Parsed phylogenetic tree object.

    Returns:
        pd.DataFrame: MARKED_dv with additional columns:
            - IN_CN_1: True if SNV overlaps CN1 region in any sample
            - IN_CN_0: True if SNV overlaps CN0 region in any sample
            - CN_1_SUBLINES: Samples with CN1 at this position
            - CN_0_SUBLINES: Samples with CN0 at this position
            - IN_CN_1_0: True if in either CN0 or CN1 (union for regenotyping)
            - CN_1_0_SUBLINES: Combined list of CN0 and CN1 samples
            - CN_1_0_MRCA_TERMINALS: Leaf nodes under CN_1_0 MRCA
            - CN_OTHER_HAPLO_BOOL: True if SNV and CNA share samples (parsimony)

    Notes:
        - CN0 = copy number 0 (homozygous deletion, no copies remain)
        - CN1 = copy number 1 (heterozygous deletion, one copy remains)
        - Both CN0 and CN1 provide regenotyping evidence (combined into CN_1_0)
        - Parsimony check is same as for Severus - blocks regenotyping if
          sample has both the SNV and the CNA event
    """

    # Create interval trees for CNA data
    # Three sets of trees: all CNA, CN1-only, and CN0-only
    # Trees are keyed by "{sample}-{chromosome}" for sample-specific lookups
    wakhan_cna_trees_per_chromosome = {}          # All copy number states
    wakhan_cna_1_only_trees_per_chromosome = {}   # CN1 (heterozygous deletion)
    wakhan_cna_0_only_trees_per_chromosome = {}   # CN0 (homozygous deletion)

    # Use provided chromosomes or fall back to autosomes for backward compatibility
    chrom_list = chromosomes if chromosomes else th_utils.autosomes

    for sub in sample_list:
        for chrom in chrom_list:
            wakhan_cna_trees_per_chromosome[f"{sub}-{chrom}"] = th_utils.it.IntervalTree()
            wakhan_cna_1_only_trees_per_chromosome[f"{sub}-{chrom}"] = th_utils.it.IntervalTree()
            wakhan_cna_0_only_trees_per_chromosome[f"{sub}-{chrom}"] = th_utils.it.IntervalTree()

    for subline in sample_list:
        bed_path = th_utils.find_cna_file(cna_path, subline)
        wk_copy_num = th_utils.read_cna_file(bed_path)
        wk_copy_num['chr'] = wk_copy_num['chr'].astype(str)
        wk_copy_num = th_utils.keep_rows_by_values(wk_copy_num, 'chr', chrom_list)
        wk_copy_num['copynumber_state'] = wk_copy_num['copynumber_state'].astype(int)

        for index, row in wk_copy_num.iterrows():
            wakhan_cna_trees_per_chromosome[f"{subline}-{row['chr']}"].addi(
                int(row['start']),
                int(row['end']) + 1,
                (str(subline), row['copynumber_state'], row['coverage'],
                 row['confidence'], row['svs_breakpoints_ids'])
            )

            if row['copynumber_state'] == 1:
                wakhan_cna_1_only_trees_per_chromosome[f"{subline}-{row['chr']}"].addi(
                    int(row['start']),
                    int(row['end']) + 1,
                    (str(subline), row['copynumber_state'], row['coverage'],
                     row['confidence'], row['svs_breakpoints_ids'])
                )

            if row['copynumber_state'] == 0:
                wakhan_cna_0_only_trees_per_chromosome[f"{subline}-{row['chr']}"].addi(
                    int(row['start']),
                    int(row['end']) + 1,
                    (str(subline), row['copynumber_state'], row['coverage'],
                     row['confidence'], row['svs_breakpoints_ids'])
                )

    # Helper functions
    def in_copy_num_of_1_helper(row):
        for subline in sample_list:
            if wakhan_cna_1_only_trees_per_chromosome[f"{subline}-{row['CHROM']}"].overlaps_point(int(row['POS'])):
                return True
        return False

    def in_copy_num_of_0_helper(row):
        for subline in sample_list:
            if wakhan_cna_0_only_trees_per_chromosome[f"{subline}-{row['CHROM']}"].overlaps_point(int(row['POS'])):
                return True
        return False

    def copy_num_1_sublines_helper(row):
        output = []
        for subline in sample_list:
            if wakhan_cna_1_only_trees_per_chromosome[f"{subline}-{row['CHROM']}"].overlaps_point(int(row['POS'])):
                output.append(str(subline))
        return output

    def copy_num_0_sublines_helper(row):
        output = []
        for subline in sample_list:
            if wakhan_cna_0_only_trees_per_chromosome[f"{subline}-{row['CHROM']}"].overlaps_point(int(row['POS'])):
                output.append(str(subline))
        return output

    # Add CNA metadata to MARKED_dv
    df_merged_copy_wakhan = MARKED_dv.copy(deep=True)

    df_merged_copy_wakhan['IN_CN_1'] = df_merged_copy_wakhan.apply(in_copy_num_of_1_helper, axis=1)
    df_merged_copy_wakhan['IN_CN_0'] = df_merged_copy_wakhan.apply(in_copy_num_of_0_helper, axis=1)
    df_merged_copy_wakhan['CN_1_SUBLINES'] = df_merged_copy_wakhan.apply(copy_num_1_sublines_helper, axis=1)
    df_merged_copy_wakhan['CN_0_SUBLINES'] = df_merged_copy_wakhan.apply(copy_num_0_sublines_helper, axis=1)

    df_merged_copy_wakhan['CN_1_MRCA'] = df_merged_copy_wakhan.apply(
        lambda row: th_utils.common_ancestor_helper(row, "CN_1_SUBLINES", input_tree=imported_tree, tree_metadata=tree_metadata)
            if row['IN_CN_1'] else float("nan"),
        axis=1
    )
    df_merged_copy_wakhan['CN_0_MRCA'] = df_merged_copy_wakhan.apply(
        lambda row: th_utils.common_ancestor_helper(row, "CN_0_SUBLINES", input_tree=imported_tree, tree_metadata=tree_metadata)
            if row['IN_CN_0'] else float("nan"),
        axis=1
    )

    df_merged_copy_wakhan['CN_1_MRCA_TERMINALS'] = df_merged_copy_wakhan.apply(
        lambda row: imported_tree.search_nodes(name=row['CN_1_MRCA'])[0].get_leaf_names()
            if row['IN_CN_1'] else float("nan"),
        axis=1
    )
    df_merged_copy_wakhan['CN_0_MRCA_TERMINALS'] = df_merged_copy_wakhan.apply(
        lambda row: imported_tree.search_nodes(name=row['CN_0_MRCA'])[0].get_leaf_names()
            if row['IN_CN_0'] else float("nan"),
        axis=1
    )

    # CN 1 and 0 Combined for final union regenotyping
    df_merged_copy_wakhan['IN_CN_1_0'] = df_merged_copy_wakhan.apply(
        lambda row: row['IN_CN_1'] or row['IN_CN_0'], axis=1
    )
    df_merged_copy_wakhan['CN_1_0_SUBLINES'] = df_merged_copy_wakhan.apply(
        lambda row: row['CN_1_SUBLINES'] + row['CN_0_SUBLINES'], axis=1
    )
    df_merged_copy_wakhan['CN_1_0_MRCA'] = df_merged_copy_wakhan.apply(
        lambda row: th_utils.common_ancestor_helper(row, "CN_1_0_SUBLINES", input_tree=imported_tree, tree_metadata=tree_metadata)
            if row['IN_CN_1_0'] else float("nan"),
        axis=1
    )
    df_merged_copy_wakhan['CN_1_0_MRCA_TERMINALS'] = df_merged_copy_wakhan.apply(
        lambda row: imported_tree.search_nodes(name=row['CN_1_0_MRCA'])[0].get_leaf_names()
            if row['IN_CN_1_0'] else float("nan"),
        axis=1
    )

    # Apply parsimony assumption for CNAs
    df_merged_copy_wakhan['CN_OTHER_HAPLO_BOOL'] = df_merged_copy_wakhan.apply(
        lambda row: len(set(row['DV_SUBLINES']).intersection(set(row['CN_1_0_SUBLINES']))) > 0
            if row['IN_CN_1_0'] else float("nan"),
        axis=1
    )
    df_merged_copy_wakhan['CN_OTHER_HAPLO_BOOL'] = df_merged_copy_wakhan['CN_OTHER_HAPLO_BOOL'].astype('boolean')

    # Merge CNA columns back to MARKED_dv
    cna_cols = ['CHROM', 'POS', 'IN_CN_1_0', 'CN_1_0_SUBLINES', 'CN_1_0_MRCA_TERMINALS',
                'IN_CN_1', 'IN_CN_0', 'CN_0_SUBLINES', 'CN_1_SUBLINES',
                'CN_0_MRCA_TERMINALS', 'CN_1_MRCA_TERMINALS', 'CN_OTHER_HAPLO_BOOL']
    cna_to_merge = df_merged_copy_wakhan[cna_cols].copy(deep=True)

    MARKED_dv = MARKED_dv.merge(cna_to_merge, on=['CHROM', 'POS'], how='left')

    return MARKED_dv


# =============================================================================
# REGENOTYPING LOGIC
# =============================================================================

def apply_regenotyping_logic(MARKED_dv, regenotype_mode, min_support_req, sample_list, imported_tree, tree_metadata):
    """
    Apply regenotyping logic based on the selected mode.

    This is the dispatcher function that calculates original DV metadata and then
    calls the appropriate mode-specific regenotyping function. All modes implement
    "dramatic shift" prevention to avoid radical changes in variant placement.

    Args:
        MARKED_dv (pd.DataFrame): SNV DataFrame with SV/CNA metadata added.
        regenotype_mode (str): One of 'severus', 'cna', or 'both'.
        min_support_req (dict): Lookup table mapping clade_size -> minimum samples.
        sample_list (list): List of sample names.
        imported_tree (ete3.Tree): Parsed phylogenetic tree object.

    Returns:
        pd.DataFrame: MARKED_dv with regenotyping columns added, plus FINAL_MRCA
            column containing the final placement decision.

    Notes:
        All modes calculate these columns:
        - DV_SUBLINES: Original samples where variant was called
        - DV_MRCA: Original placement based on called samples
        - DV_MINIMUM_SUPPORT_MET: Whether original placement meets threshold
        - REGENO_*: Mode-specific regenotyping columns
        - FINAL_MRCA: Final placement (regenotyped or original)
    """

    # Calculate original DV metadata
    MARKED_dv['DV_SUBLINES'] = MARKED_dv.apply(
        lambda row: [col for col in sample_list
                    if not col.startswith("CN") and pd.notna(row[col])],
        axis=1
    )
    MARKED_dv['DV_MRCA'] = MARKED_dv.apply(
        lambda row: th_utils.common_ancestor_helper(row, "DV_SUBLINES", input_tree=imported_tree, tree_metadata=tree_metadata),
        axis=1
    )
    MARKED_dv['DV_MRCA_TERMINALS'] = MARKED_dv.apply(
        lambda row: imported_tree.search_nodes(name=row['DV_MRCA'])[0].get_leaf_names(),
        axis=1
    )
    MARKED_dv['DV_SUBLINE_COUNT'] = MARKED_dv.apply(
        lambda row: len(row['DV_SUBLINES']), axis=1
    )
    MARKED_dv['DV_MINIMUM_SUPPORT_MET'] = MARKED_dv.apply(
        lambda row: len(row['DV_SUBLINES']) >= min_support_req[len(row['DV_MRCA_TERMINALS'])],
        axis=1
    )

    # Apply regenotyping based on mode
    if regenotype_mode == 'severus':
        MARKED_dv = apply_severus_only_regenotyping(MARKED_dv, min_support_req, sample_list, imported_tree, tree_metadata)
    elif regenotype_mode == 'cna':
        MARKED_dv = apply_cna_only_regenotyping(MARKED_dv, min_support_req, sample_list, imported_tree, tree_metadata)
    elif regenotype_mode == 'both':
        MARKED_dv = apply_union_regenotyping(MARKED_dv, min_support_req, sample_list, imported_tree, tree_metadata)

    # Perform category scoring and get final placements
    MARKED_dv = perform_category_scoring(MARKED_dv, regenotype_mode, min_support_req, sample_list, imported_tree, tree_metadata)

    return MARKED_dv


def apply_severus_only_regenotyping(MARKED_dv, min_support_req, sample_list, imported_tree, tree_metadata):
    """
    Apply regenotyping using only Severus SV deletion data.

    For each SNV that overlaps a Severus deletion, this function:
    1. Creates a union of DV_SUBLINES and SEVERUS_SUBLINES (expanded sublines)
    2. Calculates the new MRCA based on expanded sublines
    3. Checks for "dramatic shift" (>2x expansion) and blocks if detected
    4. Evaluates support threshold for the expanded placement

    The dramatic shift check prevents variants from being moved to radically
    different tree positions due to potentially spurious SV evidence.

    Args:
        MARKED_dv (pd.DataFrame): SNV DataFrame with Severus metadata.
        min_support_req (dict): Clade size -> minimum samples lookup.
        sample_list (list): List of sample names.
        imported_tree (ete3.Tree): Parsed phylogenetic tree object.

    Returns:
        pd.DataFrame: MARKED_dv with columns:
            - REGENO_COMBINED_SUBLINES_SEVERUS_ONLY: Union of DV and Severus sublines
            - REGENO_MRCA_SEVERUS_ONLY: MRCA of combined sublines
            - REGENO_MRCA_SEVERUS_ONLY_DIFFERS: Whether MRCA changed from original
            - DRAMATIC_SHIFT_SEVERUS: True if combined > 2x original sublines
            - REGENO_SEVERUS_ONLY_MINIMUM_SUPPORT_MET: Support check (False if dramatic)

    Notes:
        Regenotyping only applies when SEVERUS_OTHER_HAPLO_BOOL=False (parsimony
        assumption is satisfied - SNV and deletion are on same haplotype).
    """

    # Union of DV and Severus deletion sublines
    MARKED_dv['REGENO_COMBINED_SUBLINES_SEVERUS_ONLY'] = MARKED_dv.apply(
        lambda row: list(set(row['SEVERUS_SUBLINES']).union(set(row['DV_SUBLINES'])))
            if (row.get('IN_SEVERUS_DELETION') == True and row.get('SEVERUS_OTHER_HAPLO_BOOL') == False)
            else float("nan"),
        axis=1
    )

    # Common ancestor
    MARKED_dv['REGENO_MRCA_SEVERUS_ONLY'] = MARKED_dv.apply(
        lambda row: th_utils.common_ancestor_helper(row, 'REGENO_COMBINED_SUBLINES_SEVERUS_ONLY', input_tree=imported_tree, tree_metadata=tree_metadata)
            if (row.get('IN_SEVERUS_DELETION') and row.get('SEVERUS_OTHER_HAPLO_BOOL') == False)
            else float("nan"),
        axis=1
    )

    # Check if MRCA differs
    MARKED_dv['REGENO_MRCA_SEVERUS_ONLY_DIFFERS'] = MARKED_dv.apply(
        lambda row: row['REGENO_MRCA_SEVERUS_ONLY'] != row['DV_MRCA']
            if (row.get('IN_SEVERUS_DELETION') and row.get('SEVERUS_OTHER_HAPLO_BOOL') == False)
            else float("nan"),
        axis=1
    )
    MARKED_dv['REGENO_MRCA_SEVERUS_ONLY_DIFFERS'] = MARKED_dv['REGENO_MRCA_SEVERUS_ONLY_DIFFERS'].astype('boolean')

    # Terminal nodes
    MARKED_dv['REGENO_MRCA_SEVERUS_ONLY_TERMINALS'] = MARKED_dv.apply(
        lambda row: imported_tree.search_nodes(name=row['REGENO_MRCA_SEVERUS_ONLY'])[0].get_leaf_names()
            if row.get('IN_SEVERUS_DELETION') and row.get('SEVERUS_OTHER_HAPLO_BOOL') == False
            else float("nan"),
        axis=1
    )

    # Check for dramatic shift - prevents radical tree position changes
    # A "dramatic shift" occurs when regenotyped sublines > 2x original sublines
    # This threshold prevents a variant originally in 3 samples from being placed
    # in a clade of 10+ samples just because of SV evidence
    def tree_shift_helper_severus(row):
        if isinstance(row['REGENO_COMBINED_SUBLINES_SEVERUS_ONLY'], list):
            if len(row['REGENO_COMBINED_SUBLINES_SEVERUS_ONLY']) > (len(row['DV_SUBLINES']) * 2):
                return True
        return False

    MARKED_dv['DRAMATIC_SHIFT_SEVERUS'] = MARKED_dv.apply(
        lambda row: tree_shift_helper_severus(row)
            if row.get('IN_SEVERUS_DELETION') and row.get('SEVERUS_OTHER_HAPLO_BOOL') == False
            else float("nan"),
        axis=1
    )
    MARKED_dv['DRAMATIC_SHIFT_SEVERUS'] = MARKED_dv['DRAMATIC_SHIFT_SEVERUS'].astype('boolean')

    # Support threshold check - only passes if:
    # 1. Variant is in a Severus deletion (IN_SEVERUS_DELETION=True)
    # 2. Parsimony check passes (SEVERUS_OTHER_HAPLO_BOOL=False)
    # 3. No dramatic shift detected (DRAMATIC_SHIFT_SEVERUS=False)
    # 4. Combined sublines meet threshold for the new clade size
    MARKED_dv['REGENO_SEVERUS_ONLY_MINIMUM_SUPPORT_MET'] = MARKED_dv.apply(
        lambda row: (len(row['REGENO_COMBINED_SUBLINES_SEVERUS_ONLY']) >=
                    min_support_req[len(row['REGENO_MRCA_SEVERUS_ONLY_TERMINALS'])])
            if (row.get('IN_SEVERUS_DELETION') and row.get('SEVERUS_OTHER_HAPLO_BOOL') == False
                and row.get('DRAMATIC_SHIFT_SEVERUS') == False)
            else float("nan"),
        axis=1
    )
    MARKED_dv['REGENO_SEVERUS_ONLY_MINIMUM_SUPPORT_MET'] = MARKED_dv['REGENO_SEVERUS_ONLY_MINIMUM_SUPPORT_MET'].astype('boolean')

    return MARKED_dv


def apply_cna_only_regenotyping(MARKED_dv, min_support_req, sample_list, imported_tree, tree_metadata):
    """
    Apply regenotyping using only CNA (copy number alteration) data.

    For each SNV that overlaps a CN0 or CN1 region, this function:
    1. Creates a union of DV_SUBLINES and CN_1_0_SUBLINES (expanded sublines)
    2. Calculates the new MRCA based on expanded sublines
    3. Checks for "dramatic shift" (>2x expansion) and blocks if detected
    4. Evaluates support threshold for the expanded placement

    Args:
        MARKED_dv (pd.DataFrame): SNV DataFrame with CNA metadata.
        min_support_req (dict): Clade size -> minimum samples lookup.
        sample_list (list): List of sample names.
        imported_tree (ete3.Tree): Parsed phylogenetic tree object.

    Returns:
        pd.DataFrame: MARKED_dv with columns:
            - REGENO_COMBINED_SUBLINES_CN_1_0_ONLY: Union of DV and CNA sublines
            - REGENO_MRCA_CN_1_0_ONLY: MRCA of combined sublines
            - REGENO_MRCA_CN_1_0_ONLY_DIFFERS: Whether MRCA changed from original
            - DRAMATIC_SHIFT_CNA: True if combined > 2x original sublines
            - REGENO_CN_1_0_ONLY_MINIMUM_SUPPORT_MET: Support check (False if dramatic)

    Notes:
        Regenotyping only applies when CN_OTHER_HAPLO_BOOL=False (parsimony
        assumption is satisfied - SNV and CNA are on same haplotype).
    """

    # Union of DV and CNA sublines
    MARKED_dv['REGENO_COMBINED_SUBLINES_CN_1_0_ONLY'] = MARKED_dv.apply(
        lambda row: list(set(row['CN_1_0_SUBLINES']).union(set(row['DV_SUBLINES'])))
            if row.get('IN_CN_1_0') and row.get('CN_OTHER_HAPLO_BOOL') == False
            else float("nan"),
        axis=1
    )

    # Common ancestor
    MARKED_dv['REGENO_MRCA_CN_1_0_ONLY'] = MARKED_dv.apply(
        lambda row: th_utils.common_ancestor_helper(row, 'REGENO_COMBINED_SUBLINES_CN_1_0_ONLY', input_tree=imported_tree, tree_metadata=tree_metadata)
            if row.get('IN_CN_1_0') and row.get('CN_OTHER_HAPLO_BOOL') == False
            else float("nan"),
        axis=1
    )

    # Check if MRCA differs
    MARKED_dv['REGENO_MRCA_CN_1_0_ONLY_DIFFERS'] = MARKED_dv.apply(
        lambda row: row['REGENO_MRCA_CN_1_0_ONLY'] != row['DV_MRCA']
            if row.get('IN_CN_1_0') and row.get('CN_OTHER_HAPLO_BOOL') == False
            else float("nan"),
        axis=1
    )
    MARKED_dv['REGENO_MRCA_CN_1_0_ONLY_DIFFERS'] = MARKED_dv['REGENO_MRCA_CN_1_0_ONLY_DIFFERS'].astype('boolean')

    # Terminal nodes
    MARKED_dv['REGENO_MRCA_CN_1_0_ONLY_TERMINALS'] = MARKED_dv.apply(
        lambda row: imported_tree.search_nodes(name=row['REGENO_MRCA_CN_1_0_ONLY'])[0].get_leaf_names()
            if row.get('IN_CN_1_0') and row.get('CN_OTHER_HAPLO_BOOL') == False
            else float("nan"),
        axis=1
    )

    # Check for dramatic shift - same 2x threshold as Severus mode
    def tree_shift_helper_cna(row):
        if isinstance(row['REGENO_COMBINED_SUBLINES_CN_1_0_ONLY'], list):
            if len(row['REGENO_COMBINED_SUBLINES_CN_1_0_ONLY']) > (len(row['DV_SUBLINES']) * 2):
                return True
        return False

    MARKED_dv['DRAMATIC_SHIFT_CNA'] = MARKED_dv.apply(
        lambda row: tree_shift_helper_cna(row)
            if row.get('IN_CN_1_0') and row.get('CN_OTHER_HAPLO_BOOL') == False
            else float("nan"),
        axis=1
    )
    MARKED_dv['DRAMATIC_SHIFT_CNA'] = MARKED_dv['DRAMATIC_SHIFT_CNA'].astype('boolean')

    # Support threshold check with dramatic shift prevention
    MARKED_dv['REGENO_CN_1_0_ONLY_MINIMUM_SUPPORT_MET'] = MARKED_dv.apply(
        lambda row: (len(row['REGENO_COMBINED_SUBLINES_CN_1_0_ONLY']) >=
                    min_support_req[len(row['REGENO_MRCA_CN_1_0_ONLY_TERMINALS'])])
            if (row.get('IN_CN_1_0') and row.get('CN_OTHER_HAPLO_BOOL') == False
                and row.get('DRAMATIC_SHIFT_CNA') == False)
            else float("nan"),
        axis=1
    )
    MARKED_dv['REGENO_CN_1_0_ONLY_MINIMUM_SUPPORT_MET'] = MARKED_dv['REGENO_CN_1_0_ONLY_MINIMUM_SUPPORT_MET'].astype('boolean')

    return MARKED_dv


def apply_union_regenotyping(MARKED_dv, min_support_req, sample_list, imported_tree, tree_metadata):
    """
    Apply regenotyping using both Severus SV and CNA data combined.

    This is the default regenotyping mode ('both'). It creates
    a union of all evidence sources: DV_SUBLINES + SEVERUS_SUBLINES + CN_1_0_SUBLINES.

    For each SNV that overlaps either a Severus deletion OR a CNA region:
    1. Creates a three-way union of all sublines
    2. Calculates the new MRCA based on combined evidence
    3. Checks for "dramatic shift" (>2x expansion) and blocks if detected
    4. Evaluates support threshold for the expanded placement

    Args:
        MARKED_dv (pd.DataFrame): SNV DataFrame with both Severus and CNA metadata.
        min_support_req (dict): Clade size -> minimum samples lookup.
        sample_list (list): List of sample names.
        imported_tree (ete3.Tree): Parsed phylogenetic tree object.

    Returns:
        pd.DataFrame: MARKED_dv with columns:
            - REGENO_COMBINED_SUBLINES_UNION: Union of DV, Severus, and CNA sublines
            - REGENO_MRCA_UNION: MRCA of combined sublines
            - REGENO_MRCA_UNION_DIFFERS: Whether MRCA changed from original
            - DRAMATIC_SHIFT: True if combined > 2x original sublines
            - REGENO_UNION_MINIMUM_SUPPORT_MET: Support check (False if dramatic)

    Notes:
        - Regenotyping applies if EITHER SV or CNA parsimony check passes
        - The union can include evidence from both sources simultaneously
        - This provides the most comprehensive regenotyping but is also the most
          aggressive in terms of potential tree position changes
    """

    def combined_sublines_union_helper(row):
        dv_set = set(row['DV_SUBLINES'])
        CN_set = set()
        SEV_set = set()
        if row.get('IN_CN_1_0'):
            CN_set = set(row['CN_1_0_SUBLINES'])
        if row.get('IN_SEVERUS_DELETION'):
            SEV_set = set(row['SEVERUS_SUBLINES'])
        if SEV_set or CN_set:
            return list(dv_set.union(CN_set).union(SEV_set))
        else:
            return float("nan")

    # Union of all sublines
    MARKED_dv['REGENO_COMBINED_SUBLINES_UNION'] = MARKED_dv.apply(combined_sublines_union_helper, axis=1)

    # Common ancestor
    MARKED_dv['REGENO_MRCA_UNION'] = MARKED_dv.apply(
        lambda row: th_utils.common_ancestor_helper(row, 'REGENO_COMBINED_SUBLINES_UNION', input_tree=imported_tree, tree_metadata=tree_metadata)
            if ((row.get('IN_CN_1_0') and row.get('CN_OTHER_HAPLO_BOOL') == False) or
                (row.get('IN_SEVERUS_DELETION') and row.get('SEVERUS_OTHER_HAPLO_BOOL') == False))
            else float("nan"),
        axis=1
    )

    # Check if MRCA differs
    MARKED_dv['REGENO_MRCA_UNION_DIFFERS'] = MARKED_dv.apply(
        lambda row: row['REGENO_MRCA_UNION'] != row['DV_MRCA']
            if ((row.get('IN_CN_1_0') and row.get('CN_OTHER_HAPLO_BOOL') == False) or
                (row.get('IN_SEVERUS_DELETION') and row.get('SEVERUS_OTHER_HAPLO_BOOL') == False))
            else float("nan"),
        axis=1
    )
    MARKED_dv['REGENO_MRCA_UNION_DIFFERS'] = MARKED_dv['REGENO_MRCA_UNION_DIFFERS'].astype('boolean')

    # Terminal nodes
    MARKED_dv['REGENO_MRCA_UNION_TERMINALS'] = MARKED_dv.apply(
        lambda row: imported_tree.search_nodes(name=row['REGENO_MRCA_UNION'])[0].get_leaf_names()
            if ((row.get('IN_CN_1_0') and row.get('CN_OTHER_HAPLO_BOOL') == False) or
                (row.get('IN_SEVERUS_DELETION') and row.get('SEVERUS_OTHER_HAPLO_BOOL') == False))
            else float("nan"),
        axis=1
    )

    # Check for dramatic shift
    def tree_shift_helper_union(row):
        regeno_sublines = row['REGENO_COMBINED_SUBLINES_UNION']
        # Check if it's a list (not NaN)
        if isinstance(regeno_sublines, list):
            if len(regeno_sublines) > (len(row['DV_SUBLINES']) * 2):
                return True
        return False

    MARKED_dv['DRAMATIC_SHIFT'] = MARKED_dv.apply(
        lambda row: tree_shift_helper_union(row)
            if ((row.get('IN_CN_1_0') and row.get('CN_OTHER_HAPLO_BOOL') == False) or
                (row.get('IN_SEVERUS_DELETION') and row.get('SEVERUS_OTHER_HAPLO_BOOL') == False))
            else float("nan"),
        axis=1
    )
    MARKED_dv['DRAMATIC_SHIFT'] = MARKED_dv['DRAMATIC_SHIFT'].astype('boolean')

    # Support threshold with dramatic shift prevention
    MARKED_dv['REGENO_UNION_MINIMUM_SUPPORT_MET'] = MARKED_dv.apply(
        lambda row: (len(row['REGENO_COMBINED_SUBLINES_UNION']) >=
                    min_support_req[len(row['REGENO_MRCA_UNION_TERMINALS'])])
            if ((row.get('IN_CN_1_0') and row.get('CN_OTHER_HAPLO_BOOL') == False) or
                (row.get('IN_SEVERUS_DELETION') and row.get('SEVERUS_OTHER_HAPLO_BOOL') == False))
                and row.get('DRAMATIC_SHIFT') == False
            else float("nan"),
        axis=1
    )
    MARKED_dv['REGENO_UNION_MINIMUM_SUPPORT_MET'] = MARKED_dv['REGENO_UNION_MINIMUM_SUPPORT_MET'].astype('boolean')

    return MARKED_dv


# =============================================================================
# CATEGORY SCORING
# =============================================================================

def perform_category_scoring(MARKED_dv, regenotype_mode, min_support_req, sample_list, imported_tree, tree_metadata):
    """
    Perform category scoring and determine final variant placements.

    This function classifies each variant into one of 8 categories based on three
    boolean dimensions:
        1. Original support: Did the variant meet threshold with original DV_SUBLINES?
        2. Regenotyped support: Does it meet threshold with expanded sublines?
        3. MRCA change: Did regenotyping change the placement node?

    The 8 categories (T=True, F=False):
        Cat 1: Original=T, Regeno=T, MRCA_change=F  -> PLACED (original placement confirmed)
        Cat 2: Original=T, Regeno=T, MRCA_change=T  -> PLACED (regenotyped to different node)
        Cat 3: Original=F, Regeno=T, MRCA_change=F  -> PLACED (rescued by regenotyping)
        Cat 4: Original=F, Regeno=T, MRCA_change=T  -> PLACED (rescued and moved)
        Cat 5: Original=F, Regeno=F, MRCA_change=F  -> NOT PLACED (no support either way)
        Cat 6: Original=F, Regeno=F, MRCA_change=T  -> NOT PLACED (no support either way)
        Cat 7: Original=T, Regeno=F, MRCA_change=F  -> PLACED (use original placement)
        Cat 8: Original=T, Regeno=F, MRCA_change=T  -> PLACED (use original placement)

    Only categories 1-4 use the regenotyped placement. Categories 5-8 either fail
    placement entirely (5,6) or use the original DV placement (7,8).

    Args:
        MARKED_dv (pd.DataFrame): SNV DataFrame with all regenotyping columns.
        regenotype_mode (str): One of 'severus', 'cna', or 'both'.
        min_support_req (dict): Clade size -> minimum samples lookup.
        sample_list (list): List of sample names.
        imported_tree (ete3.Tree): Parsed phylogenetic tree object.

    Returns:
        pd.DataFrame: Combined DataFrame with FINAL_MRCA column containing the
            final placement decision, and sample columns marked with "REGENO"
            where regenotyping added evidence.

    Notes:
        - Categories 1-4 are "regenotyped" placements (use REGENO_MRCA)
        - Categories 7-8 had original support but regenotyping failed/changed nothing
        - Samples added via regenotyping are marked "REGENO" in their column
    """

    # Determine which columns to use based on mode
    if regenotype_mode == 'severus':
        regeno_col = 'REGENO_SEVERUS_ONLY_MINIMUM_SUPPORT_MET'
        mrca_col = 'REGENO_MRCA_SEVERUS_ONLY_DIFFERS'
        combined_col = 'REGENO_COMBINED_SUBLINES_SEVERUS_ONLY'
        final_mrca_col = 'REGENO_MRCA_SEVERUS_ONLY'
    elif regenotype_mode == 'cna':
        regeno_col = 'REGENO_CN_1_0_ONLY_MINIMUM_SUPPORT_MET'
        mrca_col = 'REGENO_MRCA_CN_1_0_ONLY_DIFFERS'
        combined_col = 'REGENO_COMBINED_SUBLINES_CN_1_0_ONLY'
        final_mrca_col = 'REGENO_MRCA_CN_1_0_ONLY'
    else:  # both
        regeno_col = 'REGENO_UNION_MINIMUM_SUPPORT_MET'
        mrca_col = 'REGENO_MRCA_UNION_DIFFERS'
        combined_col = 'REGENO_COMBINED_SUBLINES_UNION'
        final_mrca_col = 'REGENO_MRCA_UNION'

    # Create masks for category scoring
    original_support_TRUE = (MARKED_dv['DV_MINIMUM_SUPPORT_MET'] == True)
    original_support_FALSE = (MARKED_dv['DV_MINIMUM_SUPPORT_MET'] == False)
    regeno_support_TRUE = (MARKED_dv[regeno_col] == True)
    regeno_support_FALSE = (MARKED_dv[regeno_col] == False)
    mrca_change_TRUE = (MARKED_dv[mrca_col] == True)
    mrca_change_FALSE = (MARKED_dv[mrca_col] == False)

    # Category scoring (8 categories based on 3 boolean dimensions)
    # Categories 1-4: Regenotyping support = TRUE (these use REGENO_MRCA for placement)
    # Categories 5-8: Regenotyping support = FALSE (use original or fail)
    cat1_dv = MARKED_dv[(original_support_TRUE) & (regeno_support_TRUE) & (mrca_change_FALSE)]   # Confirmed
    cat2_dv = MARKED_dv[(original_support_TRUE) & (regeno_support_TRUE) & (mrca_change_TRUE)]    # Moved by regeno
    cat3_dv = MARKED_dv[(original_support_FALSE) & (regeno_support_TRUE) & (mrca_change_FALSE)]  # Rescued
    cat4_dv = MARKED_dv[(original_support_FALSE) & (regeno_support_TRUE) & (mrca_change_TRUE)]   # Rescued + moved
    cat5_dv = MARKED_dv[(original_support_FALSE) & (regeno_support_FALSE) & (mrca_change_FALSE)] # Failed both
    cat6_dv = MARKED_dv[(original_support_FALSE) & (regeno_support_FALSE) & (mrca_change_TRUE)]  # Failed both
    cat7_dv = MARKED_dv[(original_support_TRUE) & (regeno_support_FALSE) & (mrca_change_FALSE)]  # Use original
    cat8_dv = MARKED_dv[(original_support_TRUE) & (regeno_support_FALSE) & (mrca_change_TRUE)]   # Use original

    cat_dvs = [cat1_dv, cat2_dv, cat3_dv, cat4_dv, cat5_dv, cat6_dv, cat7_dv, cat8_dv]

    print(f"\nCategory breakdown for {regenotype_mode} mode:")
    for i, cat_df in enumerate(cat_dvs):
        print(f"Category {i+1}: {len(cat_df)} variants")

    # Add REGENO sublines column for regenotyped variants
    for i in range(len(cat_dvs)):
        if not cat_dvs[i].empty:
            # Check if combined_col has any list values (not all NaN)
            has_valid_data = cat_dvs[i][combined_col].apply(lambda x: isinstance(x, list)).any()
            if has_valid_data:
                cat_dvs[i] = cat_dvs[i].copy()
                cat_dvs[i]['ADDED_REGENO_SUBLINES'] = cat_dvs[i].apply(
                    lambda row: set(row[combined_col]) - set(row['DV_SUBLINES'])
                        if isinstance(row[combined_col], list) else set(),
                    axis=1
                )

    # Combine passing categories (1-4) - these are the regenotyped placements
    FINAL_VAR_DF = pd.concat([cat_dvs[0], cat_dvs[1], cat_dvs[2], cat_dvs[3]], ignore_index=True)

    # Mark samples that were added via regenotyping (not originally called)
    # These samples will have "REGENO" in their column instead of VCF data
    for index, row in FINAL_VAR_DF.iterrows():
        if 'ADDED_REGENO_SUBLINES' in row and isinstance(row['ADDED_REGENO_SUBLINES'], set):
            for subline in row['ADDED_REGENO_SUBLINES']:
                if subline in FINAL_VAR_DF.columns:
                    FINAL_VAR_DF.at[index, subline] = "REGENO"

    # Prepare columns for final dataframe
    sample_cols = sample_list.copy()
    columns_to_keep = ['KEY', 'CHROM', 'POS', 'REF', 'ALT', final_mrca_col] + sample_cols

    # Ensure all columns exist
    for col in columns_to_keep:
        if col not in FINAL_VAR_DF.columns:
            FINAL_VAR_DF[col] = float('nan')

    FINAL_VAR_DF = FINAL_VAR_DF[columns_to_keep]

    # Add back non-regenotyped variants
    dv_merged_copy = MARKED_dv.copy(deep=True)
    dv_merged_copy = dv_merged_copy[~dv_merged_copy['KEY'].isin(FINAL_VAR_DF['KEY'])]

    FINAL_VAR_DF.index = FINAL_VAR_DF['KEY']
    FINAL_VAR_DF = FINAL_VAR_DF.drop(columns=['KEY'])
    dv_merged_copy.index = dv_merged_copy['KEY']
    dv_merged_copy = dv_merged_copy.drop(columns=['KEY'])

    dv_merged_copy = pd.concat([dv_merged_copy, FINAL_VAR_DF], axis=0)

    # Recalculate DV metadata for combined dataframe
    dv_merged_copy['DV_SUBLINES'] = dv_merged_copy.apply(
        lambda row: [col for col in sample_list
                    if not col.startswith("CN") and pd.notna(row[col])],
        axis=1
    )
    dv_merged_copy['DV_MRCA'] = dv_merged_copy.apply(
        lambda row: th_utils.common_ancestor_helper(row, "DV_SUBLINES",
                                                     input_tree=imported_tree, tree_metadata=tree_metadata),
        axis=1
    )
    dv_merged_copy['DV_MRCA_TERMINALS'] = dv_merged_copy.apply(
        lambda row: imported_tree.search_nodes(name=row['DV_MRCA'])[0].get_leaf_names(),
        axis=1
    )
    dv_merged_copy['DV_MINIMUM_SUPPORT_MET'] = dv_merged_copy.apply(
        lambda row: len(row['DV_SUBLINES']) >= min_support_req[len(row['DV_MRCA_TERMINALS'])],
        axis=1
    )

    # Create FINAL_MRCA column
    def final_mrca_helper(row):
        regeno_mrca = row.get(final_mrca_col)
        if pd.notna(regeno_mrca):
            return regeno_mrca
        elif row['DV_MINIMUM_SUPPORT_MET']:
            return row['DV_MRCA']
        else:
            return float('nan')

    dv_merged_copy['FINAL_MRCA'] = dv_merged_copy.apply(final_mrca_helper, axis=1)

    return dv_merged_copy


# =============================================================================
# VCF OUTPUT
# =============================================================================

def generate_vcfs(MARKED_dv, non_terminal_paths, output_path, sample_list,
                  write_exclusive_vcfs, write_cumulative_vcfs, verbose, tree_metadata,
                  vcf_metadata=None, command_string=None, version=None):
    """
    Generate exclusive and cumulative VCF files for each tree node.

    This function creates two types of VCF output:
    1. Exclusive VCFs: Variants placed exactly at each node (not inherited)
    2. Cumulative VCFs: Variants at each node plus all descendant nodes

    Exclusive VCFs are useful for identifying variants that arose at specific
    branching points. Cumulative VCFs show all variants present in each lineage.

    Args:
        MARKED_dv (pd.DataFrame): SNV DataFrame with FINAL_MRCA column.
        non_terminal_paths (dict): Maps each internal node to list of all
            descendant nodes (used for cumulative VCFs).
        output_path (str): Base output directory.
        sample_list (list): List of sample names for VCF columns.
        write_exclusive_vcfs (bool): If True, write exclusive VCF files.
        write_cumulative_vcfs (bool): If True, write cumulative VCF files.
        verbose (bool): If True, print per-node variant counts.
        tree_metadata (TreeMetadata): Metadata about tree structure and node mappings.
        vcf_metadata (dict): Optional VCF metadata from input files containing:
            - 'header': VCF header string
            - 'format_string': FORMAT field string
            - 'format_definitions': Dict of FORMAT field definitions
        command_string (str): Optional command line string for provenance.

    Returns:
        tuple: (exclusive_dfs, cumulative_dfs) - dictionaries mapping node
            names to their respective DataFrames.

    Notes:
        - VCFs follow VCF4.2 format
        - Uses input VCF header and FORMAT when available
        - Output paths: {output_path}/snv/exclusive/ and {output_path}/snv/cumulative/
    """

    # Get variants that passed placement
    dv_merged_without_unplaced = MARKED_dv.dropna(subset=['FINAL_MRCA'])

    # Create node lists dynamically from tree metadata
    internal_nodes = [n for n in non_terminal_paths.keys() if tree_metadata.is_internal_node(n)]
    private_nodes = list(tree_metadata.private_node_to_leaf.keys())
    all_nodes = internal_nodes + private_nodes

    # Create exclusive DataFrames per node
    exclusive_dfs = {}
    for node in all_nodes:
        exclusive_dfs[node] = dv_merged_without_unplaced[
            dv_merged_without_unplaced['FINAL_MRCA'] == node
        ]

    # Determine FORMAT string and empty cell pattern
    if vcf_metadata and vcf_metadata.get('format_string'):
        format_string = vcf_metadata['format_string']
        format_definitions = vcf_metadata.get('format_definitions', {})
        empty_cell = th_utils.generate_empty_vcf_cell(format_string, format_definitions)
    else:
        # Fallback for backward compatibility
        format_string = "GT:GQ:DP:AD:VAF:PL"
        empty_cell = "./.:0:0:0,0:0:0,0,0"

    # Format for VCF
    sample_cols = sample_list.copy()

    # Check if we have preserved ID/QUAL/FILTER/INFO columns
    has_preserved_cols = all(col in dv_merged_without_unplaced.columns for col in ['ID', 'QUAL', 'FILTER', 'INFO'])

    for node in exclusive_dfs:
        df = exclusive_dfs[node].copy()

        # Keep necessary columns including preserved VCF columns if available
        if has_preserved_cols:
            df = df[['CHROM', 'POS', 'REF', 'ALT', 'ID', 'QUAL', 'FILTER', 'INFO'] + sample_cols]
            # Fill missing values in preserved columns
            df['ID'] = df['ID'].fillna(".")
            df['QUAL'] = df['QUAL'].fillna(".")
            df['FILTER'] = df['FILTER'].fillna("PASS")
            df['INFO'] = df['INFO'].fillna(".")
        else:
            # Fallback: create columns with default values
            df = df[['CHROM', 'POS', 'REF', 'ALT'] + sample_cols]
            df['ID'] = "."
            df['QUAL'] = "."
            df['FILTER'] = "PASS"
            df['INFO'] = "."

        # Set FORMAT column
        df['FORMAT'] = format_string

        # Fill empty sample cells with properly formatted empty string
        df.fillna(empty_cell, inplace=True)

        # Reorder columns to VCF standard order
        df = df[['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'] + sample_cols]

        exclusive_dfs[node] = df

    # Build output header
    if vcf_metadata and vcf_metadata.get('header'):
        # Use header from input VCFs
        output_header = vcf_metadata['header']
    else:
        # Fallback to generic VCF header
        output_header = """##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position (for use with symbolic alleles)">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Conditional genotype quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Read depth for each allele">
##FORMAT=<ID=VAF,Number=A,Type=Float,Description="Variant allele fractions.">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled genotype likelihoods rounded to the closest integer">
"""

    # Add TreeHarmonizer provenance
    version_str = version if version else "unknown"
    output_header += f"##TreeHarmonizer_version={version_str}\n"
    if command_string:
        output_header += f"##TreeHarmonizer_command={command_string}\n"

    # Write exclusive VCFs
    print("\nGenerating exclusive variant placement VCFs...")
    for node in exclusive_dfs:
        if verbose:
            print(f"SNV count for node {node}: {len(exclusive_dfs[node])}")

        if write_exclusive_vcfs:
            path_prefix = f"{output_path}/snv/exclusive"
            subprocess.run(['mkdir', '-p', path_prefix])
            th_utils.write_vcf(exclusive_dfs[node], f"{path_prefix}/{node}.vcf", output_header)

    # Generate cumulative VCFs
    print("Generating cumulative variant placement VCFs...")
    cumulative_dfs = {}

    for key, value in non_terminal_paths.items():
        merged_for_key = pd.concat([exclusive_dfs[x] for x in value], ignore_index=True)
        cumulative_dfs[key] = merged_for_key

        if verbose:
            print(f"SNV count for node {key}: {len(merged_for_key)}")

        if write_cumulative_vcfs:
            path_prefix = f"{output_path}/snv/cumulative"
            subprocess.run(['mkdir', '-p', path_prefix])
            th_utils.write_vcf(cumulative_dfs[key], f"{path_prefix}/{key}.vcf", output_header)

    return exclusive_dfs, cumulative_dfs
