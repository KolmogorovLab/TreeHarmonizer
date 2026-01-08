import argparse
import os
import sys
from sv import sv_placement_runner
from snv import snv_placement_runner
from cna import cna_placement_runner
from utils import (
    parse_chromosome_list, parse_sample_list, validate_chromosomes,
    get_chromosomes_for_species, get_chromosome_sizes,
    load_chromosome_sizes_file, normalize_chromosome_name,
    SPECIES_CHROMOSOMES
)
from tree_preprocessing import normalize_tree, validate_tree_for_placement

# Software version - update this when releasing new versions
TREEHARMONIZER_VERSION = "1.0.0"

def main():
    parser = argparse.ArgumentParser(description='TreeHarmonizer main program')

    # Input options for variant placement / paths
    parser.add_argument('--snv-path', type=str,
                        help='Path to SNV data directory')
    parser.add_argument('--no-snv-placement', action='store_true', default=False,
                        help='Do not perform SNV placement (default: False)')
    parser.add_argument('--cna-path', type=str,
                        help='Path to CNA data directory')
    parser.add_argument('--no-cna-placement', action='store_true', default=False,
                        help='Do not perform CNA placement (default: False)')
    parser.add_argument('--sv-path', type=str,
                        help='Path to SV VCF file')
    parser.add_argument('--no-sv-placement', action='store_true', default=False,
                        help='Do not perform SV placement (default: False)')
    
    # Sample information options
    # Optional argument for limiting samples to a specific list of folders within the variant data directories
    parser.add_argument('--sample-list', type=str, nargs='+',
                        help='List of sample folder names to include (space or comma-separated). '
                             'Example: --sample-list Sample_1 Sample_2 Sample_3 or --sample-list Sample_1,Sample_2,Sample_3. '
                             'Default: all samples in data directories.')
    
    # Regenotyping arguments
    parser.add_argument('--disable-regenotyping', action='store_true', default=False,
                        help='Disable regenotyping step (default: False)')
    parser.add_argument('--regenotype-with-sv-only', action='store_true', default=False,
                        help='Regenotype using SVs only (default: False)')
    parser.add_argument('--regenotype-with-cna-only', action='store_true', default=False,
                        help='Regenotype using CNAs only (default: False)')
    
    # FN Rate
    parser.add_argument('--fn-rate', type=float, default=0.15,
                        help='False negative rate (default: 0.15)')
    
    # Output options for VCF writing
    parser.add_argument('--write-exclusive-vcfs', action='store_true', default=True,
                        help='Write exclusive VCFs (default: True)')
    parser.add_argument('--no-write-exclusive-vcfs', dest='write_exclusive_vcfs', action='store_false',
                        help='Do not write exclusive VCFs')
    parser.add_argument('--write-cumulative-vcfs', action='store_true', default=True,
                        help='Write cumulative VCFs (default: True)')
    parser.add_argument('--no-write-cumulative-vcfs', dest='write_cumulative_vcfs', action='store_false',
                        help='Do not write cumulative VCFs')
    parser.add_argument('--output-path', type=str, default='./th_output/',
                        help='Base output path for all results (default: ./th_output/)')
    
    # Tree structure input
    parser.add_argument('--tree-newick', type=str, nargs=1,
                        help='Newick string or path to file containing tree structure (required)')
    
    ## Miscellaneous options
    parser.add_argument('--version', action='version', version=f'TreeHarmonizer {TREEHARMONIZER_VERSION}',
                        help='Show program version and exit')
    
    parser.add_argument('--verbose', action='store_true', default=False,
                        help='Enable verbose output (default: False)')

    # Reference species and chromosome selection
    parser.add_argument('--reference-species', type=str, required=True,
                        help='Reference species/genome assembly (REQUIRED). '
                             'Use "mm10" (mouse chr1-19,X,Y), "grch38" (human chr1-22,X,Y), '
                             'or provide a path to a UCSC-format chrom.sizes file.')

    # Chromosome selection (mutually exclusive group)
    chrom_group = parser.add_mutually_exclusive_group()
    chrom_group.add_argument('--chromosomes', type=str, nargs='+',
                             help='Chromosomes to profile (space or comma-separated). '
                                  'Example: --chromosomes 1 2 3 X or --chromosomes 1,2,3,X. '
                                  'Default: all chromosomes for the selected species.')
    chrom_group.add_argument('--exclude-chromosomes', type=str, nargs='+',
                             help='Chromosomes to exclude from profiling (space or comma-separated). '
                                  'Example: --exclude-chromosomes X Y. '
                                  'Mutually exclusive with --chromosomes.')

    # For testing purposes, hardcode some arguments here
    #args = parser.parse_args('--snv-path ./data/snv/ --fn-rate 0.1 --tree-newick "(A,B,(C,D));"'.split())
    #args = parser.parse_args('--snv-path ./data/snv/ --cna-path ./data/cna/ --sv-path ./data/sv/variants.vcf --fn-rate 0.1 --tree-newick "(A,B,(C,D));"'.split())
    args = parser.parse_args()

    print("args:", args)

    # First test - Tree newick string must be provided, otherwise nothing can be done
    if not args.tree_newick:
        parser.error('Newick string tree structure not provided. Use --tree-newick <newick_string> to specify the tree structure')

    # Process --tree-newick: can be either a file path or an inline newick string
    tree_newick_input = args.tree_newick[0]

    def looks_like_file_path(s):
        """Check if string looks like a file path rather than a newick string."""
        file_extensions = ('.nwk', '.txt', '.tree', '.newick')
        if s.lower().endswith(file_extensions):
            return True
        if s.startswith(('./', '../', '/')):
            return True
        return False

    if looks_like_file_path(tree_newick_input):
        # Treat as file path
        if os.path.isfile(tree_newick_input):
            with open(tree_newick_input, 'r') as f:
                tree_newick_str = f.read().strip()
            print(f"Loaded tree from file: {tree_newick_input}")
        else:
            parser.error(f"Tree file not found: {tree_newick_input}")
    else:
        # Treat as inline newick string
        tree_newick_str = tree_newick_input

    # Output path not set explicitly, use default
    if args.output_path == './th_output/':
        print("Output path not set, using default './th_output/'")

    # Determine which variant types will be placed
    # Default: if a path is provided, placement happens unless explicitly disabled
    run_snv_placement = bool(args.snv_path) and not args.no_snv_placement
    run_cna_placement = bool(args.cna_path) and not args.no_cna_placement
    run_sv_placement = bool(args.sv_path) and not args.no_sv_placement

    # Warn about contradictory input
    if args.snv_path and args.no_snv_placement:
        print("Warning: SNV path provided but --no-snv-placement flag is set. SNVs will NOT be placed.")
    if args.cna_path and args.no_cna_placement:
        print("Warning: CNA path provided but --no-cna-placement flag is set. CNAs will NOT be placed.")
    if args.sv_path and args.no_sv_placement:
        print("Warning: SV path provided but --no-sv-placement flag is set. SVs will NOT be placed.")

    # Ensure at least one variant type will be placed
    if not run_snv_placement and not run_sv_placement and not run_cna_placement:
        parser.error('No variant types will be placed! At least one of --snv-path, --cna-path, or --sv-path must be provided without its corresponding --no-placement flag.')

    # Regenotyping booleans set here.
    # Regenotyping is enabled by default, but requires:
    # 1. SNV placement to be enabled (regenotyping is part of SNV placement)
    # 2. At least one of CNA or SV paths to be provided

    # Check for conflicting regenotyping flags
    if args.disable_regenotyping and (args.regenotype_with_sv_only or args.regenotype_with_cna_only):
        parser.error('Cannot use --disable-regenotyping with --regenotype-with-sv-only or --regenotype-with-cna-only!')
    if args.regenotype_with_sv_only and args.regenotype_with_cna_only:
        parser.error('Cannot set both regenotype with SV only and regenotype with CNA only flags!')

    do_regenotyping = not args.disable_regenotyping
    if do_regenotyping and not run_snv_placement:
        do_regenotyping = False
        print("SNV placement disabled, regenotyping disabled.")
    if do_regenotyping and not args.cna_path and not args.sv_path:
        do_regenotyping = False
        print("No CNA or SV path provided, regenotyping disabled.")
    if args.regenotype_with_sv_only and not args.sv_path:
        parser.error('Regenotype with SV only flag set but no SV path provided!')
    if args.regenotype_with_cna_only and not args.cna_path:
        parser.error('Regenotype with CNA only flag set but no CNA path provided!')
    if (args.regenotype_with_sv_only or args.regenotype_with_cna_only) and not run_snv_placement:
        parser.error('Regenotype flags require SNV placement to be enabled!')

    # Print summary of all arguments applied
    print("\nTreeHarmonizer starting with the following settings:")

    print(f"SNV path: {args.snv_path}")
    print(f"CNA path: {args.cna_path}")
    print(f"SV path: {args.sv_path}")

    # Sample list handling (parse to handle both space and comma separation)
    sample_list = parse_sample_list(args.sample_list)
    if sample_list:
        print(f"Limiting samples to the following list: \n{sample_list}")
    else:
        print("Including all samples in the provided data directories.")

    # Process and validate reference species / chromosome sizes
    species_input = args.reference_species

    # Detect if input is a file path or predefined species
    if os.path.isfile(species_input):
        # Load custom chromosome sizes from file
        try:
            custom_chrom_sizes = load_chromosome_sizes_file(species_input)
        except ValueError as e:
            parser.error(str(e))
        all_species_chroms = list(custom_chrom_sizes.keys())
        species = 'custom'  # For display purposes
        is_custom_species = True
        print(f"Loaded custom chromosome sizes from: {species_input}")
        print(f"  Found {len(all_species_chroms)} chromosomes: {all_species_chroms}")
    elif species_input in SPECIES_CHROMOSOMES:
        # Use predefined species
        species = species_input
        all_species_chroms = get_chromosomes_for_species(species, include_sex=True)
        is_custom_species = False
    else:
        parser.error(
            f"Invalid --reference-species: '{species_input}' is not a valid species "
            f"(mm10, grch38) and is not a valid file path."
        )

    # Process chromosome selection arguments
    if args.chromosomes:
        # Parse the chromosome list (handles both space and comma separation)
        chromosomes_to_use = parse_chromosome_list(args.chromosomes)
        # Normalize user-provided chromosome names (strip 'chr' prefix if present)
        chromosomes_to_use = [normalize_chromosome_name(c) for c in chromosomes_to_use]

        if is_custom_species:
            # Validate against chromosomes in custom file
            invalid = set(chromosomes_to_use) - set(all_species_chroms)
            if invalid:
                parser.error(
                    f"Invalid chromosomes: {sorted(invalid)}. "
                    f"Valid chromosomes from file: {sorted(all_species_chroms)}"
                )
        else:
            # Validate against predefined species chromosomes
            try:
                chromosomes_to_use = validate_chromosomes(chromosomes_to_use, species)
            except ValueError as e:
                parser.error(str(e))
    elif args.exclude_chromosomes:
        # Start with all chromosomes and remove excluded ones
        excluded = parse_chromosome_list(args.exclude_chromosomes)
        # Normalize excluded chromosome names
        excluded = [normalize_chromosome_name(c) for c in excluded]

        if is_custom_species:
            # Validate excluded chromosomes against custom file
            invalid = set(excluded) - set(all_species_chroms)
            if invalid:
                parser.error(
                    f"Invalid chromosomes to exclude: {sorted(invalid)}. "
                    f"Valid chromosomes from file: {sorted(all_species_chroms)}"
                )
        else:
            try:
                excluded = validate_chromosomes(excluded, species)
            except ValueError as e:
                parser.error(str(e))

        chromosomes_to_use = [c for c in all_species_chroms if c not in excluded]
        if not chromosomes_to_use:
            parser.error('All chromosomes have been excluded! At least one chromosome must remain.')
    else:
        # Default: all chromosomes including X and Y (or all from custom file)
        chromosomes_to_use = all_species_chroms

    # Build chromosome sizes for downstream use
    if is_custom_species:
        if args.chromosomes or args.exclude_chromosomes:
            # Filter to selected chromosomes
            chrom_sizes = {k: v for k, v in custom_chrom_sizes.items() if k in chromosomes_to_use}
        else:
            chrom_sizes = custom_chrom_sizes.copy()
    else:
        chrom_sizes = get_chromosome_sizes(species, chromosomes_to_use)

    # Validate chromosome sizes coverage for CNA placement
    if run_cna_placement:
        missing_sizes = set(chromosomes_to_use) - set(chrom_sizes.keys())
        if missing_sizes:
            parser.error(
                f"CNA placement requires chromosome sizes for all analyzed chromosomes. "
                f"Missing sizes for: {sorted(missing_sizes)}. "
                f"Please ensure your chromosome sizes file includes these chromosomes."
            )

    print(f"Reference species: {species}" + (f" (from {species_input})" if is_custom_species else ""))
    print(f"Chromosomes to profile ({len(chromosomes_to_use)}): {chromosomes_to_use}")

    print(f"FN rate: {args.fn_rate}")
    print("Tree Newick:", tree_newick_str)

    # Validate and normalize the input tree
    try:
        validate_tree_for_placement(tree_newick_str)
    except ValueError as e:
        parser.error(f"Invalid tree: {e}")

    # Preprocess tree to ensure required structure (extra root, private nodes)
    normalized_newick, tree_metadata = normalize_tree(
        tree_newick_str,
        verbose=args.verbose
    )
    
    print(f"Run SNV placement: {run_snv_placement}")
    print(f"Run CNA placement: {run_cna_placement}")
    print(f"Run SV placement: {run_sv_placement}")
    print(f"Do regenotyping: {do_regenotyping}")
    if do_regenotyping:
        if args.regenotype_with_sv_only:
            print("  Regenotyping with SVs only.")
        elif args.regenotype_with_cna_only:
            print("  Regenotyping with CNAs only.")
        else:
            # Determine what's actually available for regenotyping
            has_cna = args.cna_path is not None
            has_sv = args.sv_path is not None
            if has_cna and has_sv:
                print("  Regenotyping with both CNAs and SVs.")
            elif has_sv:
                print("  Regenotyping with SVs only. (Only SV data provided.)")
            elif has_cna:
                print("  Regenotyping with CNAs only. (Only CNA data provided.)")
   
    print(f"Write exclusive VCFs: {args.write_exclusive_vcfs}")
    print(f"Write cumulative VCFs: {args.write_cumulative_vcfs}")
    print(f"Output path: {args.output_path}")


    ## Initialize and run TreeHarmonizer placement processes with the right parameters
    print("\nInitializing TreeHarmonizer placement process...")

    # Build command string for VCF provenance tracking
    command_string = ' '.join(sys.argv)

    # Run SNV placement if enabled
    if run_snv_placement:
        print("\n=== Running SNV Placement ===")

        # Determine regenotyping mode based on flags
        if do_regenotyping:
            if args.regenotype_with_sv_only:
                regenotype_mode = 'severus'
            elif args.regenotype_with_cna_only:
                regenotype_mode = 'cna'
            else:
                regenotype_mode = 'both'
        else:
            regenotype_mode = None

        snv_results = snv_placement_runner(
            snv_path=args.snv_path,
            tree_newick=normalized_newick,
            tree_metadata=tree_metadata,
            output_path=args.output_path,
            fn_rate=args.fn_rate,
            do_regenotyping=do_regenotyping,
            sv_path=args.sv_path,
            cna_path=args.cna_path,
            regenotype_mode=regenotype_mode,
            write_exclusive_vcfs=args.write_exclusive_vcfs,
            write_cumulative_vcfs=args.write_cumulative_vcfs,
            sample_list=sample_list,
            chromosomes=chromosomes_to_use,
            verbose=args.verbose,
            command_string=command_string,
            version=TREEHARMONIZER_VERSION
        )
        print("SNV placement completed.")

    # Run SV placement if enabled
    if run_sv_placement:
        print("\n=== Running SV Placement ===")
        sv_results = sv_placement_runner(
            sv_path=args.sv_path,
            tree_newick=normalized_newick,
            tree_metadata=tree_metadata,
            output_path=args.output_path,
            fn_rate=args.fn_rate,
            write_exclusive_vcfs=args.write_exclusive_vcfs,
            write_cumulative_vcfs=args.write_cumulative_vcfs,
            chromosomes=chromosomes_to_use,
            all_species_chroms=all_species_chroms,
            sample_list=sample_list,
            verbose=args.verbose,
            command_string=command_string,
            version=TREEHARMONIZER_VERSION
        )
        print("SV placement completed.")

    # Run CNA placement if enabled
    if run_cna_placement:
        print("\n=== Running CNA Placement ===")
        cna_results = cna_placement_runner(
            cna_path=args.cna_path,
            tree_newick=normalized_newick,
            tree_metadata=tree_metadata,
            output_path=args.output_path,
            fn_rate=args.fn_rate,
            sample_list=sample_list,
            chromosomes=chromosomes_to_use,
            chrom_sizes=chrom_sizes,
            verbose=args.verbose
        )
        print("CNA placement completed.")



if __name__ == "__main__":
    main()
