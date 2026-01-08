"""
TreeHarmonizer CNA Placement Module

This module handles placement of copy number alterations (CNAs) called by Wakhan
onto the phylogenetic tree.

Workflow Overview:
    1. Load Wakhan BED data and phylogenetic tree
    2. Create interval trees for amplifications (CN > 2) and losses (CN < 2)
    3. Generate intersections of CNA ranges across samples for each tree node
    4. Apply support threshold based on clade size and false negative rate
    5. Perform exclusive top-down placement (remove parent ranges from children)
    6. Calculate and output amplification/loss percentages per node and subline

Key Concepts:
    - Wakhan is a CNA caller that outputs per-sample BED files with copy number states
    - Amplifications: copy number > 2, Losses: copy number < 2
    - Support threshold: Similar to SNV/SV placement, based on fn_rate
    - Exclusive placement: Parent node ranges are removed from child nodes

Output:
    - TSV files with amplification/loss percentages per node (by chromosome + total)
    - TSV files with amplification/loss percentages per subline (by chromosome + total)
    - TSV files with average amplification/loss across all sublines
"""

import pandas as pd
import copy
import os
import utils as th_utils


# Mouse chromosome sizes (mm10/GRCm38)
MOUSE_CHROM_SIZES = {
    '1': 195471971, '2': 182113224, '3': 160039680, '4': 156508116,
    '5': 151834684, '6': 149736546, '7': 145441459, '8': 129401213,
    '9': 124595110, '10': 130694993, '11': 122082543, '12': 120129022,
    '13': 120421639, '14': 124902244, '15': 104043685, '16': 98207768,
    '17': 94987271, '18': 90702639, '19': 61431566
}

TOTAL_GENOME_SIZE = sum(MOUSE_CHROM_SIZES.values())


# =============================================================================
# DATA LOADING FUNCTIONS
# =============================================================================

def discover_samples(cna_path):
    """
    Discover sample names from the CNA data directory structure.

    Supports two directory structures:
    1. Simplified: cna_path/[sample_name]/[sample_name].bed
    2. Legacy (Wakhan): cna_path/[sample_name]/bed_output/[sample_name]_copynumbers_segments.bed

    Args:
        cna_path: Path to CNA data directory

    Returns:
        list: Sample names found in the directory
    """
    if not os.path.isabs(cna_path):
        cna_path = os.path.abspath(cna_path)

    sample_list = []
    for item in os.listdir(cna_path):
        item_path = os.path.join(cna_path, item)
        if os.path.isdir(item_path) and not item.startswith('_'):
            # Check new simplified path first
            bed_file = os.path.join(item_path, f'{item}.bed')
            if os.path.isfile(bed_file):
                sample_list.append(item)
                continue

            # Fall back to legacy Wakhan path
            bed_file = os.path.join(item_path, 'bed_output', f'{item}_copynumbers_segments.bed')
            if os.path.isfile(bed_file):
                sample_list.append(item)

    if len(sample_list) == 0:
        raise ValueError(f"No valid CNA sample directories found in: {cna_path}")

    return sorted(sample_list)


def load_cna_data(cna_path, sample_list, autosomes):
    """
    Load CNA data from BED files into interval trees.

    Supports both Wakhan BED format (7 columns) and generic CNA format (4 columns).

    Creates three sets of interval trees per sample-chromosome:
    - All CNA states
    - CN=1 only (heterozygous loss)
    - CN=0 only (homozygous loss)

    Also creates amplification-only and loss-only interval trees.

    Args:
        cna_path: Path to CNA data directory
        sample_list: List of sample names to load
        autosomes: List of autosome chromosome names

    Returns:
        tuple: (amp_only_dicts, loss_only_dicts, wakhan_cna_trees_per_chromosome)
    """
    if not os.path.isabs(cna_path):
        cna_path = os.path.abspath(cna_path)

    # Initialize interval trees
    wakhan_cna_trees_per_chromosome = {}
    wakhan_cna_trees_per_chromosome_amp_only = {}
    wakhan_cna_trees_per_chromosome_loss_only = {}

    for sub in sample_list:
        for chrom in autosomes:
            key = f"{sub}-{chrom}"
            wakhan_cna_trees_per_chromosome[key] = th_utils.it.IntervalTree()
            wakhan_cna_trees_per_chromosome_amp_only[key] = th_utils.it.IntervalTree()
            wakhan_cna_trees_per_chromosome_loss_only[key] = th_utils.it.IntervalTree()

    # Load BED files for each sample
    for subline in sample_list:
        bed_path = th_utils.find_cna_file(cna_path, subline)
        wk_copy_num = th_utils.read_cna_file(bed_path)
        wk_copy_num['chr'] = wk_copy_num['chr'].astype(str)

        # Filter to autosomes only
        wk_copy_num = th_utils.keep_rows_by_values(wk_copy_num, 'chr', autosomes)
        wk_copy_num['copynumber_state'] = wk_copy_num['copynumber_state'].astype(int)

        # Add intervals to trees
        for _, row in wk_copy_num.iterrows():
            key = f"{subline}-{row['chr']}"
            interval_data = (
                subline,
                row['copynumber_state'],
                row['coverage'],
                row['confidence'],
                row['svs_breakpoints_ids']
            )
            wakhan_cna_trees_per_chromosome[key].addi(
                int(row['start']),
                int(row['end']) + 1,
                interval_data
            )

    # Separate into amplifications and losses
    for subline in sample_list:
        for chrom in autosomes:
            key = f"{subline}-{chrom}"
            for interval in wakhan_cna_trees_per_chromosome[key]:
                cn_state = interval.data[1]
                if cn_state > 2:
                    # Amplification - store only the subline name as data
                    wakhan_cna_trees_per_chromosome_amp_only[key].addi(
                        interval.begin, interval.end, interval.data[0]
                    )
                elif cn_state < 2:
                    # Loss
                    wakhan_cna_trees_per_chromosome_loss_only[key].addi(
                        interval.begin, interval.end, interval.data[0]
                    )

    # Merge overlapping/tangential intervals
    for subline in sample_list:
        for chrom in autosomes:
            key = f"{subline}-{chrom}"
            wakhan_cna_trees_per_chromosome_amp_only[key].merge_overlaps(
                strict=False, data_reducer=lambda x, y: x
            )
            wakhan_cna_trees_per_chromosome_loss_only[key].merge_overlaps(
                strict=False, data_reducer=lambda x, y: x
            )

    # Convert to list-based dictionaries for processing
    amp_only_dicts = {}
    loss_only_dicts = {}

    for subline in sample_list:
        for chrom in autosomes:
            key = f"{subline}-{chrom}"
            final_amp_list = []
            final_loss_list = []

            for interval in wakhan_cna_trees_per_chromosome_amp_only[key].items():
                final_amp_list.append([interval.begin, interval.end, set([interval.data])])

            for interval in wakhan_cna_trees_per_chromosome_loss_only[key].items():
                # Exclude blank centromere region intervals
                if interval.begin == 0 and interval.end in [3000000, 3000001, 3150000]:
                    continue
                final_loss_list.append([interval.begin, interval.end, set([interval.data])])

            amp_only_dicts[key] = final_amp_list
            loss_only_dicts[key] = final_loss_list

    return amp_only_dicts, loss_only_dicts, wakhan_cna_trees_per_chromosome


# =============================================================================
# INTERVAL PROCESSING FUNCTIONS
# =============================================================================

def make_all_partial_overlap_fragments(interval_set_a, interval_keys_a, interval_set_b, interval_keys_b):
    """
    Create fragments from partial overlaps between two interval sets.

    Handles two cases:
    - Case 1: w-------y-_-_-_-_x______z (interval_a starts before interval_b)
    - Case 2: y-----w-_-_-_-_z_____x (interval_b starts before interval_a)

    Args:
        interval_set_a: First list of intervals [start, end, set_of_sublines]
        interval_keys_a: Set of keys for interval_set_a
        interval_set_b: Second list of intervals
        interval_keys_b: Set of keys for interval_set_b

    Returns:
        tuple: (overlapping_intervals, all_intervals, all_keys)
    """
    overlapping_intervals = []
    overlapping_keys = set()

    for interval_a in interval_set_a:
        for interval_b in interval_set_b:
            # Case 1: w-------y-_-_-_-_x______z
            if (interval_a[0] < interval_b[0] and
                interval_b[0] < interval_a[1] and
                interval_a[1] < interval_b[1]):
                new_interval = [
                    interval_b[0],
                    interval_a[1],
                    set(interval_a[2]).union(set(interval_b[2]))
                ]
                overlapping_intervals.append(new_interval)
                overlapping_keys.add((
                    interval_a[0],
                    interval_b[1],
                    frozenset(new_interval[2])
                ))

            # Case 2: y-----w-_-_-_-_z_____x
            elif (interval_b[0] < interval_a[0] and
                  interval_a[0] < interval_b[1] and
                  interval_b[1] < interval_a[1]):
                new_interval = [
                    interval_a[0],
                    interval_b[1],
                    set(interval_a[2]).union(set(interval_b[2]))
                ]
                overlapping_intervals.append(new_interval)
                overlapping_keys.add((
                    interval_a[0],
                    interval_b[1],
                    frozenset(new_interval[2])
                ))

    all_intervals = interval_set_a + interval_set_b + overlapping_intervals
    all_keys = interval_keys_a.union(interval_keys_b).union(overlapping_keys)

    return overlapping_intervals, all_intervals, all_keys


def merge_superset_subset_fragments(interval_set, interval_keys):
    """
    Merge intervals that are supersets/subsets or perfect overlaps.

    Handles:
    - Perfect overlaps: ac-_-_-_-_-_bd (same start and end)
    - Subset: c----a____b-----d (interval inside another)
    - Superset: a_____c-----d_____b (interval contains another)

    Args:
        interval_set: List of intervals [start, end, set_of_sublines]
        interval_keys: Set of interval keys

    Returns:
        tuple: (merged_interval_set, merged_interval_keys)
    """
    if len(interval_set) <= 1:
        return interval_set, interval_keys

    intervals_to_be_added = []
    intervals_to_be_removed = []
    interval_keys_to_be_added = set()
    interval_keys_to_be_removed = set()

    for x in range(len(interval_set)):
        for y in range(x + 1, len(interval_set)):
            int_x = interval_set[x]
            int_y = interval_set[y]

            # Perfect overlap
            if int_x[0] == int_y[0] and int_x[1] == int_y[1]:
                if int_x[2] == int_y[2]:
                    # Same data, remove duplicate
                    intervals_to_be_removed.append(int_y)
                    interval_keys_to_be_removed.add((int_y[0], int_y[1], frozenset(int_y[2])))
                else:
                    # Different data, merge
                    new_interval = [int_x[0], int_x[1], set(int_x[2]).union(set(int_y[2]))]
                    new_key = (new_interval[0], new_interval[1], frozenset(new_interval[2]))
                    if new_key not in interval_keys_to_be_added:
                        intervals_to_be_added.append(new_interval)
                        interval_keys_to_be_added.add(new_key)
                    intervals_to_be_removed.append(int_x)
                    intervals_to_be_removed.append(int_y)
                    interval_keys_to_be_removed.add((int_x[0], int_x[1], frozenset(int_x[2])))
                    interval_keys_to_be_removed.add((int_y[0], int_y[1], frozenset(int_y[2])))

            # int_y is subset of int_x
            elif int_x[0] <= int_y[0] and int_y[1] <= int_x[1]:
                if int_x[2] == int_y[2]:
                    intervals_to_be_removed.append(int_y)
                    interval_keys_to_be_removed.add((int_y[0], int_y[1], frozenset(int_y[2])))
                else:
                    int_y[2] = set(int_x[2]).union(set(int_y[2]))

            # int_x is subset of int_y
            elif int_y[0] <= int_x[0] and int_x[1] <= int_y[1]:
                if int_x[2] == int_y[2]:
                    intervals_to_be_removed.append(int_x)
                    interval_keys_to_be_removed.add((int_x[0], int_x[1], frozenset(int_x[2])))
                else:
                    int_x[2] = set(int_x[2]).union(set(int_y[2]))

    # Remove marked intervals
    for interval in intervals_to_be_removed:
        try:
            interval_set.remove(interval)
        except ValueError:
            continue
        try:
            interval_keys.remove((interval[0], interval[1], frozenset(interval[2])))
        except KeyError:
            continue

    # Add new merged intervals
    for interval in intervals_to_be_added:
        interval_set.append(interval)
        interval_keys.add((interval[0], interval[1], frozenset(interval[2])))

    return interval_set, interval_keys


def run_procedure(dict_list, key_list, current_index, prev_part_all, prev_part_keys):
    """
    Recursive procedure for merging intervals across multiple samples.

    Args:
        dict_list: List of interval lists (one per sample)
        key_list: List of key sets (one per sample)
        current_index: Current position in the lists
        prev_part_all: Accumulated intervals from previous iterations
        prev_part_keys: Accumulated keys from previous iterations

    Returns:
        list: Final merged interval list
    """
    if len(dict_list) == 1 or current_index == len(dict_list):
        return prev_part_all

    if current_index == len(dict_list) - 1:
        # Last node
        _, current_part_all, current_part_keys = make_all_partial_overlap_fragments(
            dict_list[current_index], key_list[current_index],
            prev_part_all, prev_part_keys
        )
        current_part_reduced, _ = merge_superset_subset_fragments(
            current_part_all, current_part_keys
        )
        return current_part_reduced

    elif current_index == 0:
        # First node
        _, current_part_all, current_part_keys = make_all_partial_overlap_fragments(
            dict_list[current_index], key_list[current_index],
            dict_list[current_index + 1], key_list[current_index + 1]
        )
        current_part_reduced, current_part_reduced_keys = merge_superset_subset_fragments(
            current_part_all, current_part_keys
        )
        return run_procedure(
            dict_list, key_list, current_index + 2,
            current_part_reduced, current_part_reduced_keys
        )

    else:
        # Middle node
        _, current_part_all, current_part_keys = make_all_partial_overlap_fragments(
            dict_list[current_index], key_list[current_index],
            prev_part_all, prev_part_keys
        )
        current_part_reduced, current_part_reduced_keys = merge_superset_subset_fragments(
            current_part_all, current_part_keys
        )
        return run_procedure(
            dict_list, key_list, current_index + 1,
            current_part_reduced, current_part_reduced_keys
        )


def generate_all_intersections_per_node_per_chrom(node, chrom, non_terminal_leaves,
                                                   amp_only_dicts, loss_only_dicts,
                                                   amp_or_loss="amp"):
    """
    Generate all interval intersections for a given node and chromosome.

    Args:
        node: Tree node name
        chrom: Chromosome name
        non_terminal_leaves: Dict mapping nodes to their leaf samples
        amp_only_dicts: Amplification interval dictionaries
        loss_only_dicts: Loss interval dictionaries
        amp_or_loss: "amp" for amplifications, "loss" for losses

    Returns:
        list: Merged intervals with sample support information
    """
    leaves = non_terminal_leaves[node]
    new_dict_list = []
    new_key_list = []

    for leaf in leaves:
        key = f"{leaf}-{chrom}"
        if amp_or_loss == "amp":
            temp_dict = copy.deepcopy(amp_only_dicts.get(key, []))
        else:
            temp_dict = copy.deepcopy(loss_only_dicts.get(key, []))
        temp_dict.sort(key=lambda x: (x[0], x[1]))
        new_dict_list.append(temp_dict)

    for interval_list in new_dict_list:
        temp_keys = set()
        for interval in interval_list:
            temp_keys.add((interval[0], interval[1], frozenset(interval[2])))
        new_key_list.append(temp_keys)

    output = run_procedure(new_dict_list, new_key_list, 0, [], set())
    return output if output else []


# =============================================================================
# SUPPORT THRESHOLD FUNCTIONS
# =============================================================================

def calculate_support_thresholds(clade_sizes, fn_rate):
    """
    Calculate minimum subline support requirements based on false negative rate.

    Args:
        clade_sizes: Set of clade sizes in the tree
        fn_rate: False negative rate (e.g., 0.15 for 15%)

    Returns:
        dict: Mapping of clade size to minimum support requirement
    """
    thresholds = {}
    for clade_size in clade_sizes:
        if clade_size < 2:
            thresholds[clade_size] = 1
        elif clade_size == 2:
            thresholds[clade_size] = 2
        else:
            thresholds[clade_size] = int(clade_size * (1 - fn_rate))
    return thresholds


def reduce_intersections_to_min_clade_support(all_intersections, non_terminal_leaves,
                                               support_thresholds, autosomes, tree_metadata):
    """
    Filter intervals to only those meeting minimum support threshold.

    Args:
        all_intersections: Dict of intervals per node-chromosome
        non_terminal_leaves: Dict mapping nodes to their leaves
        support_thresholds: Dict mapping clade size to min support
        autosomes: List of autosome names
        tree_metadata: TreeMetadata object for node type checks

    Returns:
        dict: Filtered intervals meeting support threshold
    """
    reduced = copy.deepcopy(all_intersections)

    for node in non_terminal_leaves:
        if tree_metadata.is_private_node(node):
            continue
        clade_size = len(non_terminal_leaves[node])
        min_support = support_thresholds.get(clade_size, 1)

        for chrom in autosomes:
            key = f"{node}-{chrom}"
            if key not in reduced:
                continue

            to_remove = []
            for interval in reduced[key]:
                if len(interval[2]) < min_support:
                    to_remove.append(interval)

            for interval in to_remove:
                try:
                    reduced[key].remove(interval)
                except ValueError:
                    continue

    return reduced


# =============================================================================
# PLACEMENT FUNCTIONS
# =============================================================================

def merge_overlapping_intervals(interval_set):
    """
    Merge overlapping intervals into unified ranges.

    Args:
        interval_set: List of intervals [start, end, set_of_sublines]

    Returns:
        list: Merged intervals
    """
    if len(interval_set) <= 1:
        return interval_set

    for x in range(len(interval_set)):
        for y in range(x + 1, len(interval_set)):
            if x >= len(interval_set) or y >= len(interval_set):
                break
            # Check for any overlap
            if interval_set[x][0] <= interval_set[y][1] and interval_set[y][0] <= interval_set[x][1]:
                new_interval = [
                    min(interval_set[x][0], interval_set[y][0]),
                    max(interval_set[x][1], interval_set[y][1]),
                    set(interval_set[x][2]).union(set(interval_set[y][2]))
                ]
                interval_set.pop(y)
                interval_set.pop(x)
                interval_set.append(new_interval)
                return merge_overlapping_intervals(interval_set)

    return interval_set


def remove_intervals_ranges_from_cur_node(intervals_to_exclude, target_node_interval_set):
    """
    Remove parent node ranges from child node intervals.

    Used for exclusive placement - ensures child nodes don't include
    ranges already placed at ancestor nodes.

    Args:
        intervals_to_exclude: List of intervals to subtract
        target_node_interval_set: Target intervals to modify

    Returns:
        list: Modified intervals with exclusions applied
    """
    for cur_to_exclude in intervals_to_exclude:
        for cur_target in list(target_node_interval_set):
            # Check for overlap
            if cur_to_exclude[0] <= cur_target[1] and cur_to_exclude[1] >= cur_target[0]:
                # If exclusion starts after target starts
                if cur_to_exclude[0] > cur_target[0]:
                    new_interval = [cur_target[0], cur_to_exclude[0], set(cur_target[2])]
                    target_node_interval_set.append(new_interval)
                # If exclusion ends before target ends
                if cur_to_exclude[1] < cur_target[1]:
                    new_interval = [cur_to_exclude[1], cur_target[1], set(cur_target[2])]
                    target_node_interval_set.append(new_interval)
                # Remove the original interval
                if cur_target in target_node_interval_set:
                    target_node_interval_set.remove(cur_target)

    return target_node_interval_set


def exclusive_range_placement(unionized_intersections, non_terminal_paths,
                               non_terminal_leaves, amp_only_dicts, loss_only_dicts,
                               autosomes, tree_metadata, amp_or_loss="amp"):
    """
    Perform exclusive top-down placement of CNA ranges.

    Traverses tree from root to leaves, removing parent ranges from children.

    Args:
        unionized_intersections: Merged intervals per node-chromosome
        non_terminal_paths: Dict mapping nodes to their path from root
        non_terminal_leaves: Dict mapping nodes to their leaves
        amp_only_dicts: Amplification intervals per sample
        loss_only_dicts: Loss intervals per sample
        autosomes: List of autosome names
        tree_metadata: TreeMetadata object for node type and name mappings
        amp_or_loss: "amp" or "loss"

    Returns:
        dict: Exclusive intervals per node-chromosome
    """
    to_mod = copy.deepcopy(unionized_intersections)
    final_exclusive = {}

    for target_node, path_to_target in non_terminal_paths.items():
        for chrom in autosomes:
            key = f"{target_node}-{chrom}"

            # Root node - no exclusions needed
            if target_node == list(non_terminal_paths.keys())[0]:
                final_exclusive[key] = to_mod.get(key, [])
                continue

            # Get target node intervals
            if tree_metadata.is_private_node(target_node):
                # Private node - get from original sample data
                sample_name = tree_metadata.get_leaf_for_private(target_node)
                sample_key = f"{sample_name}-{chrom}"
                if amp_or_loss == "amp":
                    target_intervals = copy.deepcopy(amp_only_dicts.get(sample_key, []))
                else:
                    target_intervals = copy.deepcopy(loss_only_dicts.get(sample_key, []))
            else:
                target_intervals = copy.deepcopy(to_mod.get(key, []))

            # Collect intervals from all ancestors
            intervals_to_exclude = []
            for ancestor_node in path_to_target:
                if ancestor_node == target_node:
                    continue
                ancestor_key = f"{ancestor_node}-{chrom}"
                if ancestor_key in final_exclusive:
                    for interval in final_exclusive[ancestor_key]:
                        intervals_to_exclude.append(interval)

            # Remove ancestor intervals from target
            reduced = remove_intervals_ranges_from_cur_node(
                copy.deepcopy(intervals_to_exclude),
                copy.deepcopy(target_intervals)
            )
            final_exclusive[key] = reduced

    return final_exclusive


# =============================================================================
# OUTPUT GENERATION FUNCTIONS
# =============================================================================

def calculate_percentages_per_node(exclusive_amp, exclusive_loss, non_terminal_leaves,
                                    autosomes, chrom_sizes):
    """
    Calculate amplification and loss percentages per node.

    Args:
        exclusive_amp: Exclusive amplification ranges per node-chromosome
        exclusive_loss: Exclusive loss ranges per node-chromosome
        non_terminal_leaves: Dict mapping nodes to their leaves
        autosomes: List of autosome names
        chrom_sizes: Dict mapping chromosome to size in bp

    Returns:
        tuple: (amp_percentages_df, loss_percentages_df)
    """
    total_genome = sum(chrom_sizes.values())
    amp_data = []
    loss_data = []

    for node in non_terminal_leaves:
        amp_row = {'Node': node}
        loss_row = {'Node': node}
        total_amp = 0
        total_loss = 0

        for chrom in autosomes:
            key = f"{node}-{chrom}"
            chrom_amp = sum(i[1] - i[0] for i in exclusive_amp.get(key, []))
            chrom_loss = sum(i[1] - i[0] for i in exclusive_loss.get(key, []))

            amp_row[f'chr{chrom}_pct'] = (chrom_amp / chrom_sizes[chrom]) * 100
            amp_row[f'chr{chrom}_bp'] = chrom_amp
            loss_row[f'chr{chrom}_pct'] = (chrom_loss / chrom_sizes[chrom]) * 100
            loss_row[f'chr{chrom}_bp'] = chrom_loss

            total_amp += chrom_amp
            total_loss += chrom_loss

        amp_row['total_pct'] = (total_amp / total_genome) * 100
        amp_row['total_bp'] = total_amp
        loss_row['total_pct'] = (total_loss / total_genome) * 100
        loss_row['total_bp'] = total_loss

        amp_data.append(amp_row)
        loss_data.append(loss_row)

    return pd.DataFrame(amp_data), pd.DataFrame(loss_data)


def calculate_percentages_per_subline(amp_only_dicts, loss_only_dicts, sample_list,
                                       autosomes, chrom_sizes):
    """
    Calculate amplification and loss percentages per subline (sample).

    Args:
        amp_only_dicts: Amplification intervals per sample-chromosome
        loss_only_dicts: Loss intervals per sample-chromosome
        sample_list: List of sample names
        autosomes: List of autosome names
        chrom_sizes: Dict mapping chromosome to size in bp

    Returns:
        tuple: (amp_percentages_df, loss_percentages_df, avg_amp_pct, avg_loss_pct)
    """
    total_genome = sum(chrom_sizes.values())
    amp_data = []
    loss_data = []

    for subline in sample_list:
        amp_row = {'Subline': subline}
        loss_row = {'Subline': subline}
        total_amp = 0
        total_loss = 0

        for chrom in autosomes:
            key = f"{subline}-{chrom}"
            chrom_amp = sum(i[1] - i[0] for i in amp_only_dicts.get(key, []))
            chrom_loss = sum(i[1] - i[0] for i in loss_only_dicts.get(key, []))

            amp_row[f'chr{chrom}_pct'] = (chrom_amp / chrom_sizes[chrom]) * 100
            amp_row[f'chr{chrom}_bp'] = chrom_amp
            loss_row[f'chr{chrom}_pct'] = (chrom_loss / chrom_sizes[chrom]) * 100
            loss_row[f'chr{chrom}_bp'] = chrom_loss

            total_amp += chrom_amp
            total_loss += chrom_loss

        amp_row['total_pct'] = (total_amp / total_genome) * 100
        amp_row['total_bp'] = total_amp
        loss_row['total_pct'] = (total_loss / total_genome) * 100
        loss_row['total_bp'] = total_loss

        amp_data.append(amp_row)
        loss_data.append(loss_row)

    amp_df = pd.DataFrame(amp_data)
    loss_df = pd.DataFrame(loss_data)

    avg_amp_pct = amp_df['total_pct'].mean()
    avg_loss_pct = loss_df['total_pct'].mean()

    return amp_df, loss_df, avg_amp_pct, avg_loss_pct


def write_output_files(output_path, node_amp_df, node_loss_df,
                       subline_amp_df, subline_loss_df, avg_amp, avg_loss):
    """
    Write all output TSV files.

    Args:
        output_path: Base output directory
        node_amp_df: Node amplification percentages DataFrame
        node_loss_df: Node loss percentages DataFrame
        subline_amp_df: Subline amplification percentages DataFrame
        subline_loss_df: Subline loss percentages DataFrame
        avg_amp: Average amplification percentage across sublines
        avg_loss: Average loss percentage across sublines
    """
    cna_output_dir = os.path.join(output_path, 'cna')
    os.makedirs(cna_output_dir, exist_ok=True)

    # Write node-level outputs
    node_amp_df.to_csv(
        os.path.join(cna_output_dir, 'node_amplification_percentages.tsv'),
        sep='\t', index=False
    )
    node_loss_df.to_csv(
        os.path.join(cna_output_dir, 'node_loss_percentages.tsv'),
        sep='\t', index=False
    )

    # Write subline-level outputs
    subline_amp_df.to_csv(
        os.path.join(cna_output_dir, 'subline_amplification_percentages.tsv'),
        sep='\t', index=False
    )
    subline_loss_df.to_csv(
        os.path.join(cna_output_dir, 'subline_loss_percentages.tsv'),
        sep='\t', index=False
    )

    # Write average summaries
    avg_df = pd.DataFrame({
        'metric': ['average_amplification_pct', 'average_loss_pct'],
        'value': [avg_amp, avg_loss]
    })
    avg_df.to_csv(
        os.path.join(cna_output_dir, 'average_percentages.tsv'),
        sep='\t', index=False
    )


# =============================================================================
# MAIN ENTRY POINT
# =============================================================================

def cna_placement_runner(cna_path, tree_newick, tree_metadata=None, output_path=None, fn_rate=0.15,
                          sample_list=None, chromosomes=None, chrom_sizes=None,
                          verbose=False):
    """
    Run CNA placement on Wakhan BED data.

    This is the main entry point for CNA placement. It loads Wakhan output,
    generates interval intersections per tree node, applies support thresholds,
    performs exclusive placement, and outputs percentage summaries.

    Args:
        cna_path: Path to Wakhan CNA output directory
        tree_newick: Newick-format string representing the phylogenetic tree
        tree_metadata: TreeMetadata object with node type information. If None,
            will be auto-generated from tree_newick using normalize_tree().
        output_path: Base output directory for results
        fn_rate: False negative rate for support threshold calculation (default 0.15)
        sample_list: Optional list of sample names to process (default: auto-discover)
        chromosomes: List of chromosomes to include (default: autosomes for backward compat)
        chrom_sizes: Dict mapping chromosome names to sizes in bp (default: MOUSE_CHROM_SIZES)
        verbose: If True, print detailed progress information

    Returns:
        dict: Results dictionary containing:
            - 'exclusive_amp': Dict of exclusive amplification ranges per node
            - 'exclusive_loss': Dict of exclusive loss ranges per node
            - 'node_amp_df': DataFrame with node amplification percentages
            - 'node_loss_df': DataFrame with node loss percentages
            - 'subline_amp_df': DataFrame with subline amplification percentages
            - 'subline_loss_df': DataFrame with subline loss percentages
            - 'avg_amp_pct': Average amplification percentage
            - 'avg_loss_pct': Average loss percentage
    """
    from tree_preprocessing import normalize_tree

    # Auto-generate tree_metadata if not provided (backward compatibility)
    if tree_metadata is None:
        tree_newick, tree_metadata = normalize_tree(tree_newick, verbose=verbose)

    print("Loading CNA data and tree...")

    # Configure pandas display
    pd.set_option('display.width', 5000)
    pd.set_option('display.max_columns', None)

    # Auto-discover samples if not provided
    if sample_list is None:
        sample_list = discover_samples(cna_path)
    print(f"Processing {len(sample_list)} samples: {sample_list}")

    # Use provided chromosomes or fall back to autosomes for backward compatibility
    autosomes = chromosomes if chromosomes else th_utils.autosomes

    # Use provided chromosome sizes or fall back to mouse sizes for backward compatibility
    used_chrom_sizes = chrom_sizes if chrom_sizes else MOUSE_CHROM_SIZES

    # Load tree
    imported_tree, _, _, non_terminal_paths, _, non_terminal_leaves, _, _ = \
        th_utils.get_tree_data(tree_newick, tree_metadata)

    if verbose:
        print(f"Tree loaded with {len(non_terminal_leaves)} non-terminal nodes")

    # Load CNA data
    print("Loading CNA interval data...")
    amp_only_dicts, loss_only_dicts, _ = load_cna_data(cna_path, sample_list, autosomes)

    # Calculate clade sizes and support thresholds
    print("Calculating support thresholds...")
    clade_sizes = set()
    for node in imported_tree.traverse():
        if not node.is_leaf():
            clade_sizes.add(len(node.get_leaves()))
    support_thresholds = calculate_support_thresholds(clade_sizes, fn_rate)

    if verbose:
        print(f"Support thresholds: {support_thresholds}")

    # Generate intersections for each node
    print("Generating interval intersections per node...")
    all_amp_intersections = {}
    all_loss_intersections = {}

    for node in non_terminal_leaves:
        if tree_metadata.is_private_node(node):
            continue
        for chrom in autosomes:
            key = f"{node}-{chrom}"
            all_amp_intersections[key] = generate_all_intersections_per_node_per_chrom(
                node, chrom, non_terminal_leaves, amp_only_dicts, loss_only_dicts, "amp"
            )
            all_loss_intersections[key] = generate_all_intersections_per_node_per_chrom(
                node, chrom, non_terminal_leaves, amp_only_dicts, loss_only_dicts, "loss"
            )
        if verbose:
            print(f"Processed node: {node}")

    # Apply support thresholds
    print("Applying support thresholds...")
    reduced_amp = reduce_intersections_to_min_clade_support(
        all_amp_intersections, non_terminal_leaves, support_thresholds, autosomes, tree_metadata
    )
    reduced_loss = reduce_intersections_to_min_clade_support(
        all_loss_intersections, non_terminal_leaves, support_thresholds, autosomes, tree_metadata
    )

    # Merge overlapping intervals
    print("Merging overlapping intervals...")
    unionized_amp = {}
    unionized_loss = {}
    for node in non_terminal_leaves:
        if tree_metadata.is_private_node(node):
            continue
        for chrom in autosomes:
            key = f"{node}-{chrom}"
            unionized_amp[key] = merge_overlapping_intervals(
                copy.deepcopy(reduced_amp.get(key, []))
            )
            unionized_loss[key] = merge_overlapping_intervals(
                copy.deepcopy(reduced_loss.get(key, []))
            )

    # Exclusive placement
    print("Performing exclusive placement...")
    exclusive_amp = exclusive_range_placement(
        unionized_amp, non_terminal_paths, non_terminal_leaves,
        amp_only_dicts, loss_only_dicts, autosomes, tree_metadata, "amp"
    )
    exclusive_loss = exclusive_range_placement(
        unionized_loss, non_terminal_paths, non_terminal_leaves,
        amp_only_dicts, loss_only_dicts, autosomes, tree_metadata, "loss"
    )

    # Calculate percentages
    print("Calculating percentages...")
    node_amp_df, node_loss_df = calculate_percentages_per_node(
        exclusive_amp, exclusive_loss, non_terminal_leaves, autosomes, used_chrom_sizes
    )
    subline_amp_df, subline_loss_df, avg_amp, avg_loss = calculate_percentages_per_subline(
        amp_only_dicts, loss_only_dicts, sample_list, autosomes, used_chrom_sizes
    )

    # Write output files
    print("Writing output files...")
    write_output_files(
        output_path, node_amp_df, node_loss_df,
        subline_amp_df, subline_loss_df, avg_amp, avg_loss
    )

    print(f"CNA placement complete. Output written to: {os.path.join(output_path, 'cna')}")

    return {
        'exclusive_amp': exclusive_amp,
        'exclusive_loss': exclusive_loss,
        'node_amp_df': node_amp_df,
        'node_loss_df': node_loss_df,
        'subline_amp_df': subline_amp_df,
        'subline_loss_df': subline_loss_df,
        'avg_amp_pct': avg_amp,
        'avg_loss_pct': avg_loss
    }
