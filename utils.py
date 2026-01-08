import pandas as pd
import io
from functools import reduce
import os
import intervaltree as it
import platform

#all_sublines = [1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24]
#all_sublines_with_c_char = ['C' + str(x) for x in all_sublines]

# =============================================================================
# CHROMOSOME CONFIGURATION
# =============================================================================


def normalize_chromosome_name(chrom):
    """
    Normalize chromosome name by stripping 'chr' prefix if present.

    This ensures consistent chromosome naming throughout the codebase,
    regardless of whether input files use 'chr1' or '1' format.

    Args:
        chrom: Chromosome name (e.g., 'chr1', '1', 'chrX', 'X')

    Returns:
        str: Normalized chromosome name without 'chr' prefix

    Examples:
        >>> normalize_chromosome_name('chr1')
        '1'
        >>> normalize_chromosome_name('chrX')
        'X'
        >>> normalize_chromosome_name('1')
        '1'
        >>> normalize_chromosome_name('X')
        'X'
    """
    chrom_str = str(chrom)
    if chrom_str.lower().startswith('chr'):
        return chrom_str[3:]
    return chrom_str


# Reference species chromosome configurations
SPECIES_CHROMOSOMES = {
    'mm10': {
        'autosomes': [str(x) for x in range(1, 20)],  # 1-19
        'sex_chromosomes': ['X', 'Y'],
        'all_chromosomes': [str(x) for x in range(1, 20)] + ['X', 'Y'],
    },
    'grch38': {
        'autosomes': [str(x) for x in range(1, 23)],  # 1-22
        'sex_chromosomes': ['X', 'Y'],
        'all_chromosomes': [str(x) for x in range(1, 23)] + ['X', 'Y'],
    }
}

# Chromosome sizes for percentage calculations (base pairs)
CHROMOSOME_SIZES = {
    'mm10': {
        '1': 195471971, '2': 182113224, '3': 160039680, '4': 156508116,
        '5': 151834684, '6': 149736546, '7': 145441459, '8': 129401213,
        '9': 124595110, '10': 130694993, '11': 122082543, '12': 120129022,
        '13': 120421639, '14': 124902244, '15': 104043685, '16': 98207768,
        '17': 94987271, '18': 90702639, '19': 61431566,
        'X': 171031299, 'Y': 91744698
    },
    'grch38': {
        '1': 248956422, '2': 242193529, '3': 198295559, '4': 190214555,
        '5': 181538259, '6': 170805979, '7': 159345973, '8': 145138636,
        '9': 138394717, '10': 133797422, '11': 135086622, '12': 133275309,
        '13': 114364328, '14': 107043718, '15': 101991189, '16': 90338345,
        '17': 83257441, '18': 80373285, '19': 58617616, '20': 64444167,
        '21': 46709983, '22': 50818468, 'X': 156040895, 'Y': 57227415
    }
}


def get_chromosomes_for_species(species, include_sex=True):
    """
    Get the chromosome list for a species.

    Args:
        species: 'mm10' or 'grch38'
        include_sex: If True, include X and Y chromosomes

    Returns:
        list: Chromosome names as strings
    """
    if species not in SPECIES_CHROMOSOMES:
        raise ValueError(f"Unknown species: {species}. Valid options: mm10, grch38")

    if include_sex:
        return SPECIES_CHROMOSOMES[species]['all_chromosomes'].copy()
    else:
        return SPECIES_CHROMOSOMES[species]['autosomes'].copy()


def get_chromosome_sizes(species, chromosomes=None):
    """
    Get chromosome sizes for a species.

    Args:
        species: 'mm10' or 'grch38'
        chromosomes: Optional list of chromosomes to filter to

    Returns:
        dict: Mapping of chromosome name to size in bp
    """
    if species not in CHROMOSOME_SIZES:
        raise ValueError(f"Unknown species: {species}. Valid options: mm10, grch38")

    sizes = CHROMOSOME_SIZES[species]
    if chromosomes:
        return {k: v for k, v in sizes.items() if k in chromosomes}
    return sizes.copy()


def load_chromosome_sizes_file(file_path):
    """
    Load chromosome sizes from a UCSC-format chrom.sizes file.

    Expected format: tab-separated (or comma-separated), two columns:
    - Column 1: Chromosome name (e.g., 'chr1' or '1')
    - Column 2: Chromosome size in base pairs

    The 'chr' prefix is stripped if present to ensure consistent naming.

    Args:
        file_path: Path to chrom.sizes file (TSV or CSV)

    Returns:
        dict: Mapping of chromosome name (without 'chr' prefix) to size in bp

    Raises:
        ValueError: If file format is invalid or no valid entries found

    Examples:
        # File content: "chr1\\t248956422\\nchr2\\t242193529"
        >>> sizes = load_chromosome_sizes_file('hg38.chrom.sizes')
        >>> sizes['1']
        248956422
    """
    if not os.path.isfile(file_path):
        raise ValueError(f"Chromosome sizes file not found: {file_path}")

    chrom_sizes = {}
    line_number = 0

    with open(file_path, 'r') as f:
        for line in f:
            line_number += 1
            line = line.strip()

            # Skip empty lines and comments
            if not line or line.startswith('#'):
                continue

            # Detect delimiter (tab or comma)
            if '\t' in line:
                parts = line.split('\t')
            elif ',' in line:
                parts = line.split(',')
            else:
                # Try whitespace as fallback
                parts = line.split()

            if len(parts) < 2:
                raise ValueError(
                    f"Invalid format at line {line_number} in {file_path}: "
                    f"expected 2 columns (chromosome, size), got {len(parts)}"
                )

            chrom_name = parts[0].strip()
            size_str = parts[1].strip()

            # Normalize chromosome name (strip 'chr' prefix)
            chrom_name = normalize_chromosome_name(chrom_name)

            # Validate size is numeric
            try:
                size = int(size_str)
            except ValueError:
                raise ValueError(
                    f"Invalid chromosome size at line {line_number} in {file_path}: "
                    f"'{size_str}' is not a valid integer"
                )

            if size <= 0:
                raise ValueError(
                    f"Invalid chromosome size at line {line_number} in {file_path}: "
                    f"size must be positive, got {size}"
                )

            # Check for duplicate chromosomes
            if chrom_name in chrom_sizes:
                raise ValueError(
                    f"Duplicate chromosome at line {line_number} in {file_path}: "
                    f"chromosome '{chrom_name}' already defined"
                )

            chrom_sizes[chrom_name] = size

    if not chrom_sizes:
        raise ValueError(f"No valid chromosome entries found in {file_path}")

    return chrom_sizes


def parse_chromosome_list(chromosome_arg):
    """
    Parse chromosome argument that can be space-separated or comma-separated.

    Args:
        chromosome_arg: List from argparse (space-separated) or single string with commas

    Returns:
        list: Chromosome names as strings, or None if input is None
    """
    if chromosome_arg is None:
        return None

    chromosomes = []
    for item in chromosome_arg:
        # Handle comma-separated values within each argument
        if ',' in item:
            chromosomes.extend([c.strip() for c in item.split(',')])
        else:
            chromosomes.append(item.strip())

    return chromosomes


def parse_sample_list(sample_arg):
    """
    Parse sample list argument that can be space-separated or comma-separated.

    Args:
        sample_arg: List from argparse (space-separated) or single string with commas

    Returns:
        list: Sample names as strings, or None if input is None
    """
    if sample_arg is None:
        return None

    samples = []
    for item in sample_arg:
        if ',' in item:
            samples.extend([s.strip() for s in item.split(',')])
        else:
            samples.append(item.strip())

    return samples


def validate_chromosomes(chromosomes, species):
    """
    Validate that all specified chromosomes are valid for the species.

    Args:
        chromosomes: List of chromosome names
        species: 'mm10' or 'grch38'

    Returns:
        list: Validated chromosome list

    Raises:
        ValueError: If invalid chromosomes are specified
    """
    valid = set(SPECIES_CHROMOSOMES[species]['all_chromosomes'])
    invalid = set(chromosomes) - valid
    if invalid:
        raise ValueError(
            f"Invalid chromosomes for {species}: {sorted(invalid)}. "
            f"Valid chromosomes: {sorted(valid)}"
        )
    return chromosomes


# Backward compatibility - deprecated, use get_chromosomes_for_species instead
autosomes = [str(x) for x in range(1, 20)]

def parse_vcf_samples(vcf_path):
    """
    Extract sample names from VCF header line.

    Args:
        vcf_path: Path to VCF file

    Returns:
        list: Sample names (columns after FORMAT)

    Raises:
        ValueError: If no header line found or malformed VCF
    """
    with open(vcf_path, 'r') as f:
        for line in f:
            if line.startswith('#CHROM'):
                # Split on tabs and get columns
                cols = line.strip().split('\t')

                # Verify we have standard VCF columns
                if len(cols) < 9:
                    raise ValueError(f"Malformed VCF header in {vcf_path}: expected at least 9 columns, found {len(cols)}")

                # Expected VCF columns
                expected_cols = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
                if cols[:9] != expected_cols:
                    raise ValueError(f"Malformed VCF header in {vcf_path}: expected standard VCF columns {expected_cols}, found {cols[:9]}")

                # Return sample columns (everything after FORMAT)
                return cols[9:]

        # If we get here, no #CHROM line was found
        raise ValueError(f"No #CHROM header line found in VCF file: {vcf_path}")


def parse_vcf_format_definitions(vcf_path):
    """
    Parse ##FORMAT metadata lines from VCF header to extract field definitions.

    Args:
        vcf_path: Path to VCF file

    Returns:
        dict: Mapping field ID to dict with 'Number', 'Type', 'Description'

    Example return:
        {
            'GT': {'Number': '1', 'Type': 'String', 'Description': 'Genotype'},
            'AD': {'Number': 'R', 'Type': 'Integer', 'Description': 'Read depth for each allele'},
            'PL': {'Number': 'G', 'Type': 'Integer', 'Description': 'Phred-scaled...'}
        }
    """
    import re
    format_definitions = {}

    # Pattern to match ##FORMAT=<ID=...,Number=...,Type=...,Description="...">
    pattern = re.compile(
        r'##FORMAT=<ID=([^,]+),Number=([^,]+),Type=([^,]+),Description="([^"]*)">'
    )

    with open(vcf_path, 'r') as f:
        for line in f:
            if line.startswith('##FORMAT='):
                match = pattern.match(line.strip())
                if match:
                    field_id = match.group(1)
                    format_definitions[field_id] = {
                        'Number': match.group(2),
                        'Type': match.group(3),
                        'Description': match.group(4)
                    }
            elif line.startswith('#CHROM'):
                # Stop reading after header line
                break

    return format_definitions


def generate_empty_vcf_cell(format_string, format_definitions, num_alleles=2, ploidy=2):
    """
    Generate an empty/missing VCF cell value based on FORMAT string and metadata.

    Args:
        format_string: Colon-separated FORMAT field (e.g., "GT:GQ:DP:AD:VAF:PL")
        format_definitions: Dict from parse_vcf_format_definitions()
        num_alleles: Number of alleles (ref + alt), default 2 for biallelic
        ploidy: Ploidy level, default 2 for diploid

    Returns:
        str: Properly formatted empty cell (e.g., "./.:0:0:0,0:0:0,0,0")

    Number field interpretation:
        - '1': single value -> '0' for Integer/Float, '.' for String
        - 'R': ref + alt alleles -> comma-separated zeros (num_alleles values)
        - 'A': alt alleles only -> comma-separated zeros (num_alleles - 1 values)
        - 'G': genotypes -> comma-separated zeros (for diploid = 3: hom-ref, het, hom-alt)
        - '.': variable -> single '0'

    GT is special-cased to always be './.:'
    """
    fields = format_string.split(':')
    empty_values = []

    for field in fields:
        if field == 'GT':
            # Genotype field - always use ./. for missing
            empty_values.append('./.')
            continue

        # Get field definition if available
        field_def = format_definitions.get(field, {})
        number = field_def.get('Number', '1')
        field_type = field_def.get('Type', 'Integer')

        # Determine the default value based on type
        if field_type == 'String':
            default_val = '.'
        else:
            default_val = '0'

        # Determine how many values based on Number
        if number == '1' or number == '.':
            empty_values.append(default_val)
        elif number == 'R':
            # One per allele (ref + alts)
            empty_values.append(','.join([default_val] * num_alleles))
        elif number == 'A':
            # One per alt allele
            count = max(1, num_alleles - 1)
            empty_values.append(','.join([default_val] * count))
        elif number == 'G':
            # One per genotype: for diploid, this is (n+1)*n/2 where n is num_alleles
            # For biallelic diploid: 3 values (0/0, 0/1, 1/1)
            genotype_count = (num_alleles * (num_alleles + 1)) // 2
            empty_values.append(','.join([default_val] * genotype_count))
        else:
            # Try to parse as integer
            try:
                count = int(number)
                empty_values.append(','.join([default_val] * count))
            except ValueError:
                # Unknown number format, use single value
                empty_values.append(default_val)

    return ':'.join(empty_values)


def read_vcf_deepvariant(path, sample=None):
    """
    Read a DeepVariant VCF file.

    Args:
        path: Path to VCF file
        sample: Sample name to use. If None, auto-detect from VCF header.

    Returns:
        tuple: (DataFrame, vcf_header, sample_name)

    Raises:
        ValueError: If auto-detection fails or multiple samples found
    """
    # Auto-detect sample if not provided
    if sample is None:
        detected_samples = parse_vcf_samples(path)
        if len(detected_samples) == 0:
            raise ValueError(f"No sample columns found in VCF: {path}")
        if len(detected_samples) > 1:
            raise ValueError(
                f"Expected single sample in DeepVariant VCF, found {len(detected_samples)}: {detected_samples}. "
                f"File: {path}"
            )
        sample = detected_samples[0]

    with open(path, 'r') as f:
        header_list = [l for l in f if l.startswith('##')]
        vcf_header = "".join(header_list)
        f.seek(0)
        lines = [l for l in f if not l.startswith('##')]

    # Build dtype dict dynamically
    dtype_dict = {
        '#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
        'QUAL': str, 'FILTER': str, 'INFO': str, 'FORMAT': str,
        sample: str
    }

    df = pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype=dtype_dict,
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})

    # Normalize chromosome names (strip 'chr' prefix if present)
    df['CHROM'] = df['CHROM'].apply(normalize_chromosome_name)

    return df, vcf_header, sample

def read_vcf_severus(path, simple_name=False):
    """
    Read a Severus SV VCF file.

    Args:
        path: Path to VCF file
        simple_name: If False, handle .haplotagged suffix stripping

    Returns:
        tuple: (DataFrame, vcf_header, sample_list)

    Raises:
        ValueError: If no samples found
    """
    # Auto-detect sample names from VCF header
    detected_samples = parse_vcf_samples(path)

    if len(detected_samples) == 0:
        raise ValueError(f"No sample columns found in Severus VCF: {path}")

    with open(path, 'r') as f:
        header_list = [l for l in f if l.startswith('##')]
        vcf_header = "".join(header_list)
        f.seek(0)
        lines = [l for l in f if not l.startswith('##')]

    # Build dtype dict dynamically for standard VCF columns
    dtype_dict = {
        '#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
        'QUAL': str, 'FILTER': str, 'INFO': str, 'FORMAT': str
    }

    # Add sample columns to dtype dict (all samples are strings)
    for sample in detected_samples:
        dtype_dict[sample] = str

    # Handle column renaming if needed
    rename_dict = {'#CHROM': 'CHROM'}

    if not simple_name:
        # Check if samples have .haplotagged suffix and set up renaming
        for sample in detected_samples:
            if sample.endswith('.haplotagged'):
                # Map from .haplotagged name to clean name
                clean_name = sample.replace('.haplotagged', '')
                rename_dict[sample] = clean_name

    df = pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype=dtype_dict,
        sep='\t'
    ).rename(columns=rename_dict)

    # Normalize chromosome names (strip 'chr' prefix if present)
    df['CHROM'] = df['CHROM'].apply(normalize_chromosome_name)

    # Get final sample list after renaming
    if not simple_name:
        sample_list = [s.replace('.haplotagged', '') for s in detected_samples]
    else:
        sample_list = detected_samples

    return df, vcf_header, sample_list

def read_bed(path, header_input='infer'):
    df = pd.read_csv(path, sep='\t', comment='#', header=header_input)
    return df

def read_bed_updated(path):
    """
    Read a Wakhan-format BED file with 7 columns.

    Handles both header formats:
    - Headers starting with '#' (treated as comments)
    - Headers without '#' prefix (detected by checking if second column is non-numeric)

    Args:
        path: Path to the BED file

    Returns:
        pd.DataFrame: DataFrame with columns: chr, start, end, coverage,
                      copynumber_state, confidence, svs_breakpoints_ids
    """
    # Peek at file to detect header without '#'
    skip_header = False
    with open(path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            # First non-comment line - check if it's a header
            first_line = line.strip()
            break

    # Check if first line looks like a Wakhan header (text, not data)
    # Header would have: chr, start, end, coverage, copynumber_state, confidence, svs_breakpoints_ids
    parts = first_line.split('\t')
    if len(parts) >= 7:
        # If second column is not numeric, it's likely a header row
        try:
            int(parts[1])
        except ValueError:
            skip_header = True

    df = pd.read_csv(
        path, sep='\t', comment='#', header=None,
        skiprows=1 if skip_header else 0,
        names=['chr', 'start', 'end', 'coverage', 'copynumber_state', 'confidence', 'svs_breakpoints_ids']
    )
    # Normalize chromosome names (strip 'chr' prefix if present)
    df['chr'] = df['chr'].apply(normalize_chromosome_name)
    return df


# =============================================================================
# CNA FILE HANDLING
# =============================================================================

def find_cna_file(cna_path, sample_name):
    """
    Find a CNA file for a sample, checking multiple possible locations.

    Searches for CNA files in the following order:
    1. {cna_path}/{sample}/{sample}.bed (new simplified path)
    2. {cna_path}/{sample}/bed_output/{sample}_copynumbers_segments.bed (legacy Wakhan path)

    Args:
        cna_path: Base path to CNA data directory
        sample_name: Name of the sample

    Returns:
        str: Path to the found CNA file

    Raises:
        FileNotFoundError: If no CNA file is found for the sample
    """
    if not os.path.isabs(cna_path):
        cna_path = os.path.abspath(cna_path)

    # Check new simplified path first
    simplified_path = os.path.join(cna_path, sample_name, f'{sample_name}.bed')
    if os.path.isfile(simplified_path):
        return simplified_path

    # Fall back to legacy Wakhan path
    legacy_path = os.path.join(cna_path, sample_name, 'bed_output', f'{sample_name}_copynumbers_segments.bed')
    if os.path.isfile(legacy_path):
        return legacy_path

    raise FileNotFoundError(
        f"No CNA file found for sample '{sample_name}'. "
        f"Checked locations:\n"
        f"  - {simplified_path}\n"
        f"  - {legacy_path}"
    )


# Column name aliases for generic CNA format (case-insensitive matching)
_CNA_COLUMN_ALIASES = {
    'chr': ['chrom', 'chr', 'chromosome'],
    'start': ['start', 'begin', 'pos_start'],
    'end': ['end', 'stop', 'pos_end'],
    'copynumber_state': ['copy_number', 'cn', 'copynumber', 'copynumber_state']
}


def _detect_cna_format(path):
    """
    Detect whether a CNA file is in Wakhan or generic format.

    Detection logic:
    - Skip lines starting with '#' (comments)
    - Read first non-comment line
    - If line has exactly 4 columns AND contains recognized column aliases → generic
    - Otherwise → Wakhan BED format

    Args:
        path: Path to the CNA file

    Returns:
        str: 'generic' or 'wakhan'
    """
    with open(path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            first_line = line.strip()
            break

    # Try tab separator first, then comma
    if '\t' in first_line:
        parts = first_line.split('\t')
    else:
        parts = first_line.split(',')

    # Generic format must have exactly 4 columns with recognized header names
    if len(parts) == 4:
        # Check if this looks like a header row with recognized column names
        all_aliases = set()
        for aliases in _CNA_COLUMN_ALIASES.values():
            all_aliases.update(a.lower() for a in aliases)

        # Count how many parts match known column aliases
        matches = sum(1 for p in parts if p.strip().lower() in all_aliases)

        # If at least 3 of 4 columns match known aliases, treat as generic format
        if matches >= 3:
            return 'generic'

    return 'wakhan'


def _find_column(df_columns, target, aliases):
    """
    Find a column in a DataFrame using a list of possible aliases.

    Args:
        df_columns: DataFrame column names
        target: Target column name (for error message)
        aliases: List of possible column name aliases

    Returns:
        str: The found column name

    Raises:
        ValueError: If no matching column is found
    """
    df_cols_lower = {c.lower(): c for c in df_columns}
    for alias in aliases:
        if alias.lower() in df_cols_lower:
            return df_cols_lower[alias.lower()]
    raise ValueError(f"Could not find {target} column. Expected one of: {aliases}")


def read_cna_generic(path):
    """
    Read a generic CNA file with flexible column names.

    Expects a CSV or TSV file with a header row containing columns for:
    - Chromosome: chrom, chr, chromosome (case-insensitive)
    - Start: start, begin, pos_start (case-insensitive)
    - End: end, stop, pos_end (case-insensitive)
    - Copy number: copy_number, cn, copynumber (case-insensitive)

    Args:
        path: Path to the CNA file

    Returns:
        pd.DataFrame: DataFrame with standardized columns: chr, start, end, coverage,
                      copynumber_state, confidence, svs_breakpoints_ids
                      (missing Wakhan columns get default values)
    """
    # Detect separator (tab or comma)
    with open(path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            first_line = line.strip()
            break

    sep = '\t' if '\t' in first_line else ','

    # Read with header
    df = pd.read_csv(path, sep=sep, comment='#')

    # Find columns using aliases
    chr_col = _find_column(df.columns, 'chromosome', _CNA_COLUMN_ALIASES['chr'])
    start_col = _find_column(df.columns, 'start', _CNA_COLUMN_ALIASES['start'])
    end_col = _find_column(df.columns, 'end', _CNA_COLUMN_ALIASES['end'])
    cn_col = _find_column(df.columns, 'copy_number', _CNA_COLUMN_ALIASES['copynumber_state'])

    # Create standardized DataFrame
    result = pd.DataFrame({
        'chr': df[chr_col].apply(normalize_chromosome_name),
        'start': df[start_col],
        'end': df[end_col],
        'coverage': 0.0,
        'copynumber_state': df[cn_col],
        'confidence': 1.0,
        'svs_breakpoints_ids': '[]'
    })

    return result


def read_cna_file(path):
    """
    Read a CNA file, auto-detecting the format.

    Supports two formats:
    1. Wakhan BED format (7 columns, positional): chr, start, end, coverage,
       copynumber_state, confidence, svs_breakpoints_ids
    2. Generic format (4 columns, header-based): chrom, start, end, copy_number
       (with flexible column name aliases)

    Args:
        path: Path to the CNA file

    Returns:
        pd.DataFrame: DataFrame with standardized columns: chr, start, end, coverage,
                      copynumber_state, confidence, svs_breakpoints_ids
    """
    format_type = _detect_cna_format(path)

    if format_type == 'generic':
        return read_cna_generic(path)
    else:
        return read_bed_updated(path)


def write_vcf(df, path, input_header):
    with open(path, 'w') as vcf:
        vcf.write(input_header)
    df.to_csv(path, sep="\t", mode='a', index=False)

    # Determine which system we are currently on, macOS or Linux, and run the appropriate sed command to replace 'CHROM' with '#CHROM'
    curren_system = platform.system()
    if curren_system == "Linux":
        os.system("sed -i 's/CHROM/#CHROM/g' " + path)
    elif curren_system == "Darwin":  # macOS
        os.system("sed -i '' 's/CHROM/#CHROM/g' " + path)
        
    #os.system("sed -i 's/CHROM/#CHROM/g' " + path)
    # MACOS FORMAT REQUIRES '' BEFORE THE STRING, CHANGE ON LINUX SERVER
    #os.system("sed -i '' 's/CHROM/#CHROM/g' " + path)

def keep_rows_by_values(df, col, values):
    return df[df[col].isin(values)]

def generate_merged_df(caller_path, predefined_sample_list=None, chromosomes=None):
    """
    Generate a merged DataFrame from per-sample VCF files.

    Args:
        caller_path: Path to directory containing sample subdirectories
        predefined_sample_list: Optional list of sample names to include
        chromosomes: Optional list of chromosomes to filter to

    Returns:
        tuple: (caller_merged, sample_list, vcf_metadata)
            - caller_merged: Merged DataFrame with all variants
            - sample_list: List of sample names
            - vcf_metadata: Dict with 'header', 'format_string', 'format_definitions', 'format_warnings'
    """
    # If caller path is a relative path, make it absolute
    if not os.path.isabs(caller_path):
        caller_path = os.path.abspath(caller_path)

    caller_path = caller_path.strip("/")
    caller_path = "/" + caller_path

    all_caller_vcfs = []
    sample_list = []

    # VCF metadata to capture from first sample
    vcf_metadata = {
        'header': None,
        'format_string': None,
        'format_definitions': None,
        'format_warnings': []
    }

    # For preserving ID/QUAL/FILTER/INFO from first occurrence of each variant
    first_occurrence_data = None

    # If using a predefined sample list, ensure all samples are accounted for
    if predefined_sample_list is not None:
        sample_list = predefined_sample_list
        for sample in sample_list:
            if not os.path.isdir(caller_path + '/' + sample):
                raise ValueError("Sample directory not found: " + caller_path + '/' + sample)

    # If general sample structure, assume each folder inside caller path is a sample name, except those that start with "_"
    # Assume within each folder there is a VCF file names [sample_name].vcf

    else:
        # List directories that do not start with '_'
        sample_list = [
            d for d in os.listdir(caller_path)
            if os.path.isdir(os.path.join(caller_path, d)) and not d.startswith('_')
        ]
        # Remove trailing slashes and keep only directory names
        sample_list = [os.path.basename(os.path.normpath(s)) for s in sample_list]
        # Make sure VCFs exist for each sample, and that there is at least one sample
        if len(sample_list) == 0:
            print("Caller path: ", caller_path)
            raise ValueError("No sample directories found in the provided caller path.")
        for sample in sample_list:
            if not os.path.isdir(caller_path + '/' + sample):
                raise ValueError("Sample directory not found: " + caller_path + '/' + sample)
            if not os.path.isfile(caller_path + '/' + sample + '/' + sample + '.vcf'):
                raise ValueError("VCF file not found for sample: " + caller_path + '/' + sample + '/' + sample + '.vcf')

    for i, sample in enumerate(sample_list):
        vcf_path = caller_path + '/' + sample + '/' + sample + '.vcf'
        caller, vcf_header, detected_sample = read_vcf_deepvariant(
            vcf_path,
            sample=sample
        )

        # Validate that detected sample matches expected
        if detected_sample != sample:
            raise ValueError(
                f"Sample name mismatch for directory '{sample}': VCF contains '{detected_sample}'. "
                f"VCF file: {vcf_path}"
            )

        # Capture metadata from first sample
        if i == 0:
            vcf_metadata['header'] = vcf_header
            vcf_metadata['format_definitions'] = parse_vcf_format_definitions(vcf_path)
            # Get FORMAT string from first data row
            if 'FORMAT' in caller.columns and len(caller) > 0:
                vcf_metadata['format_string'] = caller['FORMAT'].iloc[0]
        else:
            # Validate FORMAT consistency across samples (warn, don't error)
            if 'FORMAT' in caller.columns and len(caller) > 0:
                current_format = caller['FORMAT'].iloc[0]
                if current_format != vcf_metadata['format_string']:
                    vcf_metadata['format_warnings'].append(
                        f"FORMAT mismatch in {sample}: expected '{vcf_metadata['format_string']}', got '{current_format}'"
                    )

        caller['CHROM'] = caller.CHROM.astype(str)
        caller['POS'] = caller.POS.astype(str)
        caller['REF'] = caller.REF.astype(str)
        caller['ALT'] = caller.ALT.astype(str)
        # Filter to specified chromosomes (defaults to autosomes for backward compatibility)
        chrom_filter = chromosomes if chromosomes else autosomes
        caller = keep_rows_by_values(caller, 'CHROM', chrom_filter)
        caller['KEY'] = caller['CHROM'].astype(str) + ":" + caller['POS'].astype(str) + ":" + caller['REF'].astype(str) + ":" + caller['ALT'].astype(str)

        # Capture ID/QUAL/FILTER/INFO for first occurrence of each variant
        if i == 0:
            first_occurrence_data = caller[['KEY', 'ID', 'QUAL', 'FILTER', 'INFO']].copy()
        else:
            # Add new KEYs not seen before
            new_keys = caller[~caller['KEY'].isin(first_occurrence_data['KEY'])]
            if len(new_keys) > 0:
                first_occurrence_data = pd.concat([
                    first_occurrence_data,
                    new_keys[['KEY', 'ID', 'QUAL', 'FILTER', 'INFO']]
                ], ignore_index=True)

        # Drop columns for the main merge (same as before, but we'll re-add them after)
        caller.drop(columns=['ID', 'POS', 'CHROM', 'FILTER', 'INFO', 'ALT', 'FORMAT', 'QUAL', 'REF'], inplace=True)
        all_caller_vcfs.append(caller)

    caller_merged = reduce(lambda left, right: pd.merge(left, right, on=['KEY'], how='outer'), all_caller_vcfs)

    # Re-add the preserved ID/QUAL/FILTER/INFO columns
    if first_occurrence_data is not None:
        caller_merged = caller_merged.merge(first_occurrence_data, on='KEY', how='left')

    return caller_merged, sample_list, vcf_metadata

def generate_severus_df(severus_path, simple_name=False, chromosomes=None, sample_list=None):
    """
    Generate a Severus SV DataFrame from VCF file.

    Args:
        severus_path: Path to Severus VCF file
        simple_name: If False, handle .haplotagged suffix stripping
        chromosomes: Optional list of chromosomes to filter to (defaults to autosomes)
        sample_list: Optional list of samples to include (defaults to all samples in VCF)

    Returns:
        tuple: (DataFrame, sample_list, vcf_header)
            - DataFrame: The Severus VCF data
            - sample_list: List of sample names
            - vcf_header: VCF header string from input file
    """
    # Filter to specified chromosomes (defaults to autosomes for backward compatibility)
    str_chrom_list_to_use = chromosomes if chromosomes else autosomes

    if simple_name:
        sev_vcf, header, detected_sample_list = read_vcf_severus(severus_path, simple_name=True)
    else:
        sev_vcf, header, detected_sample_list = read_vcf_severus(severus_path)

    sev_vcf['CHROM'] = sev_vcf.CHROM.astype(str)
    sev_vcf = keep_rows_by_values(sev_vcf, 'CHROM', str_chrom_list_to_use)

    # Filter to specified samples if provided
    if sample_list is not None:
        # Validate all requested samples exist in VCF (strict mode)
        missing_samples = [s for s in sample_list if s not in detected_sample_list]
        if missing_samples:
            raise ValueError(f"Samples not found in SV VCF: {missing_samples}. "
                           f"Available samples: {detected_sample_list}")
        # Filter DataFrame to keep only VCF fixed columns + requested sample columns
        fixed_cols = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
        sev_vcf = sev_vcf[fixed_cols + sample_list]
        return sev_vcf, sample_list, header

    return sev_vcf, detected_sample_list, header

def generate_wakhan_df(wakhan_path):
    pass

#### Tree Methods ####

from ete3 import Tree
from tree_preprocessing import normalize_tree  # noqa: F401 - re-exported for convenience


def get_tree_data(newick_tree_string, tree_metadata):
    """
    Parse tree and extract path/ancestor data structures.

    Args:
        newick_tree_string: Newick format tree string (should be normalized)
        tree_metadata: TreeMetadata object from normalize_tree()

    Returns:
        Tuple of:
        - imported_tree: The ete3 Tree object
        - non_terminals: List of non-terminal node objects
        - terminals: List of terminal (leaf) node objects
        - non_terminal_paths: Dict mapping each internal node to its path from root
        - terminal_paths: Dict mapping each leaf to its ancestral path
        - non_terminal_leaves: Dict mapping internal nodes to their descendant leaves
        - terminal_paths_o_keys: Dict mapping private nodes to their ancestral paths
        - non_terminal_paths_without_N1: Variant excluding the root's direct child
    """
    imported_tree = Tree(newick_tree_string, format=1)

    # Use metadata to get root name
    root_node_name = tree_metadata.root_name

    non_terminal_paths = {}
    terminal_paths = {}
    non_terminal_leaves = {}
    terminal_paths_o_keys = {}

    # Filter non-terminals using metadata instead of hardcoded patterns
    non_terminals = [
        node for node in imported_tree.get_descendants()
        if node.name != root_node_name and not tree_metadata.is_leaf(node.name)
    ]
    terminals = imported_tree.get_leaves()

    # Get paths to every non terminal and terminal node
    for non_terminal in non_terminals:
        path = non_terminal.get_ancestors()
        reverse_path = path[::-1]
        reverse_path.append(non_terminal)
        non_terminal_paths.update({non_terminal.name: [node.name for node in reverse_path[1:]]})

    for terminal in terminals:
        clade_list = terminal.get_ancestors()
        reverse_path = clade_list[::-1]
        terminal_paths.update({terminal.name: [node.name for node in reverse_path]})
        terminal_paths_o_keys.update({reverse_path[-1].name: [node.name for node in reverse_path[1:]]})

    # Get collection of leaves for every non terminal node
    for non_terminal in non_terminals:
        non_terminal_leaves.update({non_terminal.name: [leaf.name for leaf in non_terminal.get_leaves()]})

    # Make a version of non_terminal_paths without the root's direct child
    # (previously hardcoded as "N1")
    root_node = imported_tree.get_tree_root()
    root_children = root_node.get_children()
    removal_node = root_children[0].name if root_children else None

    non_terminal_paths_without_N1 = {}
    for key, value in non_terminal_paths.items():
        if key == removal_node:
            continue
        value_copy = value.copy()
        for element in value_copy:
            if element == removal_node:
                value_copy.remove(element)
        non_terminal_paths_without_N1.update({key: value_copy})

    return imported_tree, non_terminals, terminals, non_terminal_paths, terminal_paths, non_terminal_leaves, terminal_paths_o_keys, non_terminal_paths_without_N1


def common_ancestor_helper(row, input_col, input_tree, tree_metadata):
    """
    Find the most recent common ancestor for samples in a row.

    Args:
        row: DataFrame row containing sample names
        input_col: Column name containing list of sample names
        input_tree: The ete3 Tree object
        tree_metadata: TreeMetadata object for node name mappings

    Returns:
        Name of the most recent common ancestor node
    """
    if len(row[input_col]) == 1:
        leaf_name = row[input_col][0]
        private_name = tree_metadata.get_private_for_leaf(leaf_name)
        return input_tree.get_common_ancestor([leaf_name, private_name]).name
    else:
        return input_tree.get_common_ancestor(row[input_col]).name