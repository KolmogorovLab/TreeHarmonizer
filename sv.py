"""
TreeHarmonizer SV Placement Module

This module handles placement of structural variants (SVs) called by Severus
onto the phylogenetic tree.

Workflow Overview:
    1. Load Severus VCF data and phylogenetic tree
    2. Filter sex chromosomes (X, Y) and duplicate BND breakpoint pairs
    3. Determine which samples have each SV called (variant depth > 0)
    4. Calculate MRCA (most recent common ancestor) for each SV
    5. Apply support threshold based on clade size and false negative rate
    6. Generate exclusive and cumulative VCF files per tree node

Key Concepts:
    - Severus is an SV caller that outputs merged multi-sample VCFs
    - BND (breakend) variants have paired entries (_1 and _2); we keep only _1
    - Support threshold: Similar to SNV placement, based on fn_rate
    - Unlike SNV placement, SV placement does NOT use regenotyping

Output:
    - Exclusive VCFs: SVs placed exactly at each node
    - Cumulative VCFs: SVs at each node plus all descendant nodes
"""

import pandas as pd
import subprocess
import re
import utils as th_utils


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def normalize_breakend_chromosomes(alt_value):
    """
    Normalize chromosome references in SV breakend ALT field.

    SV breakends in VCF ALT field contain chromosome references like:
    - N[chr5:12345[  -> N[5:12345[
    - ]chrX:131935953] -> ]X:131935953]
    - [chr1:100[ -> [1:100[

    This function strips the 'chr' prefix from these references to ensure
    consistent chromosome naming throughout the codebase.

    Args:
        alt_value: ALT field value from VCF (string)

    Returns:
        str: ALT value with 'chr' prefixes removed from breakend references
    """
    # Pattern matches: [ or ] followed by 'chr' and then alphanumeric chromosome + colon
    # Captures the bracket and chromosome name, replaces with bracket + chromosome (no chr)
    return re.sub(r'([\[\]])chr([0-9XYxy]+):', r'\1\2:', str(alt_value))


def severus_called_sublines_helper(row, sample_list):
    """
    Determine which sublines have a called SV based on variant depth (DV) > 0.

    Severus VCF FORMAT field structure: GT:DR:DV:VAF:hVAF
    - GT: Genotype
    - DR: Reference read depth
    - DV: Variant read depth (index 4, 0-based)
    - VAF: Variant allele frequency
    - hVAF: Haplotype-specific VAF

    A sample is considered to have the SV if DV > 0.

    Args:
        row (pd.Series): DataFrame row containing sample columns.
        sample_list (list): List of sample column names to check.

    Returns:
        list: Sample names where this SV has variant read support (DV > 0).
    """
    output_subline_list = []
    for col in sample_list:
        # Parse FORMAT field: GT:DR:DV:VAF:hVAF (DV is at index 4)
        internal_sev_data = row[col].split(":")
        DV = int(internal_sev_data[4])
        if DV > 0:
            output_subline_list.append(col)
    return output_subline_list


def generate_exclusive_presence_absence_vcf(internal_node, input_merged_df):
    """
    Generate a DataFrame of variants exclusively placed at a given node.

    Filters the input DataFrame to return only SVs that:
    1. Meet the minimum subline support threshold
    2. Have their MRCA exactly at the specified node

    Args:
        internal_node (str): The tree node name to filter for (e.g., 'N5', 'O12').
        input_merged_df (pd.DataFrame): DataFrame with SV placement columns including
            'minimum_subline_support_threshold_met_severus' and 'severus_mrca'.

    Returns:
        pd.DataFrame: Deep copy of filtered DataFrame containing only SVs
            placed at this node with sufficient support.
    """
    df_filtered = input_merged_df[
        (input_merged_df['minimum_subline_support_threshold_met_severus'] == True) &
        (input_merged_df['severus_mrca'] == internal_node)
    ]
    return df_filtered.copy(deep=True)


# =============================================================================
# MAIN ENTRY POINT
# =============================================================================

def sv_placement_runner(sv_path, tree_newick, tree_metadata=None, output_path=None, fn_rate=0.15,
                        write_exclusive_vcfs=True, write_cumulative_vcfs=True,
                        chromosomes=None, all_species_chroms=None, sample_list=None, verbose=False,
                        command_string=None, version=None):
    """
    Run SV placement on Severus VCF data.

    This is the main entry point for SV placement. It loads Severus output,
    filters out problematic variants (breakends to excluded chromosomes, duplicate BND pairs),
    calculates MRCA for each SV, applies support thresholds, and generates
    output VCF files.

    Args:
        sv_path (str): Path to Severus merged VCF file.
        tree_newick (str): Newick-format string representing the phylogenetic tree.
        tree_metadata: TreeMetadata object with node type information. If None,
            will be auto-generated from tree_newick using normalize_tree().
        output_path (str): Base output directory for results and VCFs.
        fn_rate (float): False negative rate for support threshold calculation.
            Default 0.15 means we expect 15% of true SVs to be missed.
        write_exclusive_vcfs (bool): If True, write VCFs with SVs exclusively
            placed at each tree node.
        write_cumulative_vcfs (bool): If True, write VCFs with SVs at each
            node plus all descendant nodes.
        chromosomes (list): List of chromosomes to include in analysis.
            If None, defaults to autosomes.
        all_species_chroms (list): List of all chromosomes for the species.
            Used to determine which chromosomes are excluded for breakend filtering.
        sample_list (list): List of sample names to include in analysis.
            If None, includes all samples found in the VCF.
        verbose (bool): If True, print detailed per-node SV counts.
        command_string (str): Command string used to invoke TreeHarmonizer,
            used for provenance tracking in output VCF headers.

    Returns:
        dict: Results dictionary containing:
            - 'exclusive_dfs': Dict mapping node names to exclusive SV DataFrames
            - 'cumulative_dfs': Dict mapping node names to cumulative SV DataFrames

    Notes:
        - Writes placed/unplaced SV TSV files to {output_path}/sv/
        - Writes exclusive VCFs to {output_path}/sv/exclusive/
        - Writes cumulative VCFs to {output_path}/sv/cumulative/
        - SVs with breakends to excluded chromosomes are filtered out
        - BND variants have paired entries; only the first (_1) is kept
    """
    from tree_preprocessing import normalize_tree

    # Auto-generate tree_metadata if not provided (backward compatibility)
    if tree_metadata is None:
        tree_newick, tree_metadata = normalize_tree(tree_newick, verbose=verbose)

    print("Loading Severus data and tree...")

    # Configure pandas display options
    pd.set_option('display.width', 5000)
    pd.set_option("display.expand_frame_repr", True)
    pd.set_option("display.max_colwidth", 1000)
    pd.set_option('display.max_columns', None)

    # Load Severus data into a merged DataFrame
    severus, sample_list, sv_header = th_utils.generate_severus_df(
        severus_path=sv_path, simple_name=True, chromosomes=chromosomes, sample_list=sample_list
    )

    # Normalize chromosome references in breakend ALT field (strip 'chr' prefix)
    # This ensures breakend filtering works correctly regardless of input format
    severus['ALT'] = severus['ALT'].apply(normalize_breakend_chromosomes)

    # Load tree input via newick string and parse it into various components
    imported_tree, _, _, non_terminal_paths, _, non_terminal_leaves, _, _ = \
        th_utils.get_tree_data(tree_newick, tree_metadata)

    print("Filtering breakends to excluded chromosomes and duplicate BND pairs...")

    # Set CHROM column to type string
    severus['CHROM'] = severus['CHROM'].astype(str)

    # Determine which chromosomes are excluded (for breakend filtering)
    # Breakends pointing to excluded chromosomes must be filtered out
    if chromosomes and all_species_chroms:
        excluded_chroms = set(all_species_chroms) - set(chromosomes)
    else:
        # Backward compatibility: if no chromosome info provided, filter X and Y
        excluded_chroms = {'X', 'Y'}

    if excluded_chroms:
        # Build regex pattern to match excluded chromosomes in breakend notation
        # Breakends appear as: ]CHROM:POS] or [CHROM:POS[ or N[CHROM:POS[ etc.
        # The chromosome always appears before a colon in the ALT field
        patterns = [f"{chrom}:" for chrom in excluded_chroms]
        pattern_regex = '|'.join(patterns)

        # Find SVs with breakends to excluded chromosomes
        sev_excluded_endpoint = severus[
            severus['ALT'].str.contains(pattern_regex, regex=True, na=False)
        ].copy()

        # Also remove the paired breakend entry (change _1 to _2 in ID)
        ids_to_remove = sev_excluded_endpoint['ID'].tolist()
        id_pairs_to_remove = [x[:-1] + "2" for x in ids_to_remove if x.endswith("1")]
        all_ids_to_remove = list(set(ids_to_remove + id_pairs_to_remove))
        severus = severus[~severus['ID'].isin(all_ids_to_remove)].copy()

        if len(all_ids_to_remove) > 0:
            print(f"Filtered {len(all_ids_to_remove)} SVs with breakends to excluded chromosomes: {sorted(excluded_chroms)}")

    # Remove duplicate BND (breakend) entries
    # Severus outputs BND variants as pairs: ID_1 and ID_2 for each breakpoint
    # We only keep _1 entries to avoid double-counting the same SV
    bnd_second_breakends = severus[
        (severus['ID'].str.contains("BND") == True) &
        (severus['ID'].str.endswith("2") == True)
    ].copy()
    severus = severus[~severus['ID'].isin(bnd_second_breakends['ID'])].copy()

    sev_df = severus.copy(deep=True)

    print("Preprocessing common ancestors...")

    # Determine which sublines have each SV called
    sev_df['severus_called_sublines'] = sev_df.apply(
        lambda row: severus_called_sublines_helper(row, sample_list), axis=1
    )

    # Find the most recent common ancestor (MRCA) for each SV
    sev_df['severus_mrca'] = sev_df.apply(
        lambda row: th_utils.common_ancestor_helper(row, "severus_called_sublines", imported_tree, tree_metadata),
        axis=1
    )

    # Get terminal nodes for each MRCA
    sev_df['severus_mrca_terminals'] = sev_df.apply(
        lambda row: imported_tree.search_nodes(name=row['severus_mrca'])[0].get_leaf_names(),
        axis=1
    )

    print("Determining clade size acceptance requirements...")

    # Get all clade sizes in the tree (number of leaf nodes under each internal node)
    clade_sizes = set()
    for clade in imported_tree.traverse():
        if clade.is_leaf():
            continue
        clade_sizes.add(len(clade.get_leaves()))

    # Calculate minimum support requirements based on false negative rate
    # This determines how many samples must have an SV for it to be placed at each clade
    # Formula accounts for expected false negatives in SV calling
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

    print("Determining placed variants...")

    # Calculate support metrics for each SV
    # called_subline_count: How many samples have this SV
    # terminal_count: How many leaf nodes are under the MRCA (clade size)
    sev_df['called_subline_count'] = sev_df.apply(
        lambda row: len(row['severus_called_sublines']), axis=1
    )
    sev_df['terminal_count'] = sev_df.apply(
        lambda row: len(row['severus_mrca_terminals']), axis=1
    )

    # Check if SV meets the support threshold for its clade
    # An SV is "placed" if called_sublines >= required_support for that clade size
    sev_df['minimum_subline_support_threshold_met_severus'] = sev_df.apply(
        lambda row: row['called_subline_count'] >=
            minimum_subline_support_per_clade_size_requirement[row['terminal_count']],
        axis=1
    )

    #print("\nPlacement threshold results:")
    #print(sev_df['minimum_subline_support_threshold_met_severus'].value_counts())

    # Create output files of placed and unplaced variants
    placed_severus_df = sev_df[
        sev_df['minimum_subline_support_threshold_met_severus'] == True
    ].copy()
    unplaced_severus_df = sev_df[
        sev_df['minimum_subline_support_threshold_met_severus'] == False
    ].copy()
    placed_severus_df.to_csv(f"{output_path}/sv/placed_severus_variants.tsv", sep="\t", index=False)
    unplaced_severus_df.to_csv(f"{output_path}/sv/unplaced_severus_variants.tsv", sep="\t", index=False)
    print(f"Placed SV variants: {len(placed_severus_df)}")
    print(f"Unplaced SV variants: {len(unplaced_severus_df)}")


    # Default Severus VCF header (fallback if input header not available)
    default_severus_header = """##fileformat=VCFv4.2
##source=Severus_v1.3
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=INS,Description="Insertion">
##ALT=<ID=DUP,Description="Duplication">
##ALT=<ID=INV,Description="Reciprocal Inversion">
##ALT=<ID=BND,Description="Breakend">
##FILTER=<ID=PASS,Description="All filters passed">
##INFO=<ID=PRECISE,Number=0,Type=Flag,Description="SV with precise breakpoints coordinates and length">
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="SV with imprecise breakpoints coordinates and length">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of the SV">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the SV">
##INFO=<ID=STRANDS,Number=1,Type=String,Description="Breakpoint strandedness">
##INFO=<ID=DETAILED_TYPE,Number=1,Type=String,Description="Detailed type of the SV">
##INFO=<ID=INSLEN,Number=1,Type=Integer,Description="Length of the unmapped sequence between breakpoint">
##INFO=<ID=MAPQ,Number=1,Type=Integer,Description="Median mapping quality of supporting reads">
##INFO=<ID=MATE_ID,Number=1,Type=String,Description="MATE ID for breakends">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DR,Number=1,Type=Integer,Description="Number of reference reads">
##FORMAT=<ID=DV,Number=1,Type=Integer,Description="Number of variant reads">
##FORMAT=<ID=VAF,Number=1,Type=Float,Description="Variant allele frequency">
##FORMAT=<ID=hVAF,Number=3,Type=Float,Description="Haplotype specific variant Allele frequency (H0,H1,H2)">
"""

    # Use input header if available, otherwise use default
    if sv_header:
        output_header = sv_header
    else:
        output_header = default_severus_header

    # Add TreeHarmonizer provenance to header
    version_str = version if version else "unknown"
    output_header += f"##TreeHarmonizer_version={version_str}\n"
    if command_string:
        output_header += f"##TreeHarmonizer_command={command_string}\n"

    print("\nGenerating exclusive variant placement VCFs...")

    exclusive_dfs = {}

    # Get sample columns for reordering
    sample_columns = sample_list.copy()
    vcf_columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'] + sample_columns

    # Generate Exclusive Presence Absence VCFs for every node
    for key, value in non_terminal_leaves.items():
        framework_df_exclusive = generate_exclusive_presence_absence_vcf(key, sev_df)
        reordered_df = framework_df_exclusive[vcf_columns]
        if verbose:
            print(f"SV count for node {key}: {len(framework_df_exclusive)}")

        if write_exclusive_vcfs:
            path_prefix = f"{output_path}/sv/exclusive"
            subprocess.run(['mkdir', '-p', path_prefix])
            th_utils.write_vcf(reordered_df, f"{path_prefix}/{key}.vcf", output_header)

        exclusive_dfs[key] = framework_df_exclusive

    print("Generating cumulative variant placement VCFs...")

    cumulative_severus_dfs = {}

    for key, value in non_terminal_paths.items():
        merged_for_key = pd.concat([exclusive_dfs[x] for x in value], ignore_index=True)
        cumulative_severus_dfs[key] = merged_for_key
        if verbose:
            print(f"SV count for node {key}: {len(merged_for_key)}")

        if write_cumulative_vcfs:
            path_prefix = f"{output_path}/sv/cumulative"
            subprocess.run(['mkdir', '-p', path_prefix])
            th_utils.write_vcf(merged_for_key, f"{path_prefix}/{key}.vcf", output_header)

    return {
        'exclusive_dfs': exclusive_dfs,
        'cumulative_dfs': cumulative_severus_dfs
    }
