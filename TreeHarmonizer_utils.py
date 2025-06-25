import pandas as pd
import io
from functools import reduce
import os
import intervaltree as it

all_sublines = [1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24]
all_sublines_with_c_char = ['C' + str(x) for x in all_sublines]
chromosomes_without_13 = [str(x) for x in range(1, 20) if x != 13]
chromosomes_with_13 = [str(x) for x in range(1, 20)]

def read_vcf_deepvariant(path, sample='default'):
    with open(path, 'r') as f:
        header_list = [l for l in f if l.startswith('##')]
        vcf_header = "".join(header_list)
        f.seek(0)
        lines = [l for l in f if not l.startswith('##')]
    df = pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str, 'FORMAT': str, sample: str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})
    return df, vcf_header

def read_vcf_severus(path, simple_name=False):
    with open(path, 'r') as f:
        header_list = [l for l in f if l.startswith('##')]
        vcf_header = "".join(header_list)
        f.seek(0)
        lines = [l for l in f if not l.startswith('##')]
    if simple_name:
        df = pd.read_csv(
            io.StringIO(''.join(lines)),
            dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
                   'QUAL': str, 'FILTER': str, 'INFO': str, 'FORMAT': str, 'C1': str, 'C3': str, 'C4': str, 'C5': str, 'C6': str, 'C7': str, 'C8': str, 'C9': str, 'C10': str, 'C11': str, 'C12': str, 'C13': str, 'C14': str, 'C15': str, 'C16': str, 'C17': str, 'C18': str, 'C19': str, 'C20': str, 'C21': str, 'C22': str, 'C23': str, 'C24': str},
            sep='\t'
        ).rename(columns={'#CHROM': 'CHROM'})
    else:
        df = pd.read_csv(
            io.StringIO(''.join(lines)),
            dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
                   'QUAL': str, 'FILTER': str, 'INFO': str, 'FORMAT': str, 'C1.haplotagged': str, 'C3.haplotagged': str, 'C4.haplotagged': str, 'C5.haplotagged': str, 'C6.haplotagged': str, 'C7.haplotagged': str, 'C8.haplotagged': str, 'C9.haplotagged': str, 'C10.haplotagged': str, 'C11.haplotagged': str, 'C12.haplotagged': str, 'C13.haplotagged': str, 'C14.haplotagged': str, 'C15.haplotagged': str, 'C16.haplotagged': str, 'C17.haplotagged': str, 'C18.haplotagged': str, 'C19.haplotagged': str, 'C20.haplotagged': str, 'C21.haplotagged': str, 'C22.haplotagged': str, 'C23.haplotagged': str, 'C24.haplotagged': str},
            sep='\t'
        ).rename(columns={'#CHROM': 'CHROM', 'C1.haplotagged': 'C1', 'C3.haplotagged': 'C3', 'C4.haplotagged': 'C4', 'C5.haplotagged': 'C5', 'C6.haplotagged': 'C6', 'C7.haplotagged': 'C7', 'C8.haplotagged': 'C8', 'C9.haplotagged': 'C9', 'C10.haplotagged': 'C10', 'C11.haplotagged': 'C11', 'C12.haplotagged': 'C12', 'C13.haplotagged': 'C13', 'C14.haplotagged': 'C14', 'C15.haplotagged': 'C15', 'C16.haplotagged': 'C16', 'C17.haplotagged': 'C17', 'C18.haplotagged': 'C18', 'C19.haplotagged': 'C19', 'C20.haplotagged': 'C20', 'C21.haplotagged': 'C21', 'C22.haplotagged': 'C22', 'C23.haplotagged': 'C23', 'C24.haplotagged': 'C24'})
    return df, vcf_header

def read_bed(path, header_input='infer'):
    df = pd.read_csv(path, sep='\t', comment='#', header=header_input)
    return df

def read_bed_updated(path):
    df = pd.read_csv(path, sep='\t', comment='#', header=None, names=['chr', 'start', 'end', 'coverage', 'copynumber_state', 'confidence', 'svs_breakpoints_ids'])
    return df

def write_vcf(df, path, input_header):
    with open(path, 'w') as vcf:
        vcf.write(input_header)
    df.to_csv(path, sep="\t", mode='a', index=False)
    os.system("sed -i 's/CHROM/#CHROM/g' " + path)
    # MACOS FORMAT REQUIRES '' BEFORE THE STRING, CHANGE ON LINUX SERVER
    #os.system("sed -i '' 's/CHROM/#CHROM/g' " + path)

def keep_rows_by_values(df, col, values):
    return df[df[col].isin(values)]

def generate_merged_df(caller_name, caller_path, without_chromosome_13=False):
    caller_path.strip("/")
    caller_path = "/" + caller_path
    #print(caller_path)
    all_caller_vcfs = []
    for x in all_sublines:
        if caller_name == "dv_new":
            caller, vcf_header = read_vcf_deepvariant(caller_path + '/C' + str(x) + '/C' + str(x) + '.vcf')
        elif caller_name == "dv":
            caller, vcf_header = read_vcf_deepvariant(caller_path + '/C' + str(x) + '.vcf')
        elif caller_name == "cs":
            caller, vcf_header = read_vcf_deepvariant(caller_path + '/C' + str(x) + '/C' + str(x) + '_pass.vcf')
        caller['CHROM'] = caller.CHROM.astype(str)
        caller['POS'] = caller.POS.astype(str)
        caller['REF'] = caller.REF.astype(str)
        caller['ALT'] = caller.ALT.astype(str)
        if without_chromosome_13:
            caller = keep_rows_by_values(caller, 'CHROM', chromosomes_without_13)
        else:
            caller = keep_rows_by_values(caller, 'CHROM', chromosomes_with_13)
        caller['KEY'] = caller['CHROM'].astype(str) + ":" + caller['POS'].astype(str) + ":" + caller['REF'].astype(str) + ":" + caller['ALT'].astype(str)
        caller.drop(columns=['ID', 'POS', 'CHROM', 'FILTER', 'INFO', 'ALT', 'FORMAT', 'QUAL', 'REF'], inplace=True)

        all_caller_vcfs.append(caller)

    caller_merged = reduce(lambda  left,right: pd.merge(left,right,on=['KEY'], how='outer'), all_caller_vcfs)

    return caller_merged

def generate_severus_df(severus_path, without_chromosome_13=False, simple_name=False):
    if without_chromosome_13:
        str_chrom_list_to_use = chromosomes_without_13
    else:
        str_chrom_list_to_use = chromosomes_with_13

    if simple_name:
        sev_vcf, header = read_vcf_severus(severus_path, simple_name=True)
    else:
        sev_vcf, header = read_vcf_severus(severus_path)
    sev_vcf['CHROM'] = sev_vcf.CHROM.astype(str)
    sev_vcf = keep_rows_by_values(sev_vcf, 'CHROM', str_chrom_list_to_use)

    return sev_vcf

def generate_wakhan_df(wakhan_path):
    pass

#### Tree Methods ####

from ete3 import Tree

def get_tree_data(newick_tree_string="(((((((((((C20)O20,(C7)O7)N22,(C8)O8)N21,(C16)O16)N20,(C11)O11)N19,((C18)O18,(C15)O15)N18)N17,(C13)O13)N16,((((C21)O21,(C6)O6)N15,(C24)O24)N14,(C9)O9)N13)N12,((((C1)O1,(C22)O22)N11,(C4)O4)N10,(C5)O5)N9)N8,(((C23)O23,(((C10)O10,(C12)O12)N6,((C3)O3,(C14)O14)N7)N5)N4,((C19)O19,(C17)O17)N3)N2)N1)N0;", non_original=False):

    imported_tree = Tree(newick_tree_string, format=1)

    # If running the tree validation testing
    root_node_name = "N0"
    if non_original:
        root_node_name = imported_tree.get_tree_root().name
    #print("Root node name: ", root_node_name)
    #print(imported_tree)

    non_terminal_paths = {}
    terminal_paths = {}
    non_terminal_leaves = {}
    terminal_paths_o_keys = {}

    non_terminals = [node for node in imported_tree.get_descendants() if root_node_name not in node.name and "C" not in node.name]
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

    # Make an WITHOUT N1 version of non_terminal_paths
    if not non_original:
        removal_node = "N1"
    else:
        removal_node = "Node" + str(int(root_node_name[-2:]) - 1)
        
    non_terminal_paths_without_N1 = {}
    for key, value in non_terminal_paths.items():
        if key == removal_node:
            continue
        value_copy = value.copy()
        for element in value_copy:
            if element == removal_node:
                value_copy.remove(element)
        non_terminal_paths_without_N1.update({key: value_copy})

    #print(non_terminals)
    #print(terminals)
    #print(non_terminal_paths)
    #print(terminal_paths)
    #print(non_terminal_leaves)
    #print(terminal_paths_o_keys)
    #print(non_terminal_paths_without_N1)

    return imported_tree, non_terminals, terminals, non_terminal_paths, terminal_paths, non_terminal_leaves, terminal_paths_o_keys, non_terminal_paths_without_N1

def common_ancestor_helper(row, input_col, input_tree):
    if len(row[input_col]) == 1:
        mod_str = row[input_col][0]
        mod_str = mod_str.replace("C", "O")
        return input_tree.get_common_ancestor([row[input_col][0], mod_str]).name
    else:
        return input_tree.get_common_ancestor(row[input_col]).name