import sys
from DnDsCODes import codes
import argparse


def saveToFasta(DICT_FASTA: dict = None, output: str = None):
    with open(output, "w") as f:
        for k, s in DICT_FASTA.items():
            f.write(f">{k}\n")
            f.write(f"{s}\n")

def __get_YEARS(DICT_FASTA: dict = None):
    years = []
    for key in DICT_FASTA.keys():
        year = key.split("|")[-1].split("-")[0]
        years.append(year)
    years = set(years)
    return years

def __SPLIT_BY_YEAR(DICT_FASTA: dict = None):
    years = __get_YEARS(DICT_FASTA=DICT_FASTA)
    DICTS = {year: {} for year in years}

    for k, s in DICT_FASTA.items():
        year = k.split("|")[-1].split("-")[0]
        DICTS[year].update({k: s})
    
    DICTS_NR = { year: codes.REMOVE_REDUNDANCY(DICTS[year]) for year in years }

    return DICTS_NR

def main():

    parser = argparse.ArgumentParser(description='run_DnDsCODes', \
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--ALN_FASTA_FILE', required=True, type=str, help='', default=None)
    parser.add_argument('--TREE_FILE', required=True, type=str, help='', default=None)
    parser.add_argument('--ROOT_TREE', help='', dest='ROOT_TREE', action='store_true')
    parser.add_argument('--SPLIT_BY_YEAR', help='', dest='SPLIT_BY_YEAR', action='store_true')
    parser.add_argument('--SAVE_PHYLIP', help='', dest='SAVE_PHYLIP', action='store_true')
    
    
    args = parser.parse_args()
    
    print(args)

    DICT_FASTA = codes.READ_FASTA(FASTA_FILE=args.ALN_FASTA_FILE, NR=False)
    DICT_FASTA_WITHOUT_STOP = codes.REMOVE_STOP(DICT_FASTA)

    if args.SPLIT_BY_YEAR:
        DICTS_NR = __SPLIT_BY_YEAR(DICT_FASTA_WITHOUT_STOP)
        
        for k in DICTS_NR.keys():
            saveToFasta(DICTS_NR[k], output=f"{k}_aln_pruned_withoutStop.fasta")
            if args.SAVE_PHYLIP:
                codes.FastaToPhylip(f"{k}_aln_pruned_withoutStop.fasta", 
                        output=f"{k}_aln_pruned_withoutStop.phy")
        
            allowed_labels = list(DICTS_NR[k].keys())

            codes.PRUNE_TREE(allowed_labels=allowed_labels, 
                            TREE_FILE=args.TREE_FILE, 
                            output=f"{k}_pruned_tree.nwk")
            if args.ROOT_TREE:
                codes.ROOT_TREE(TREE_FILE=f"{k}_pruned_tree.nwk", output= f"{k}_rooted_tree.nwk")
    else:
        DICT_NR = codes.REMOVE_REDUNDANCY(DICT_FASTA_WITHOUT_STOP)
        saveToFasta(DICT_NR, output="aln_pruned_withoutStop.fasta")
        if args.SAVE_PHYLIP:
            codes.FastaToPhylip("aln_pruned_withoutStop.fasta", 
                    output="aln_pruned_withoutStop.phy")
    
        allowed_labels = list(DICT_NR.keys())

        codes.PRUNE_TREE(allowed_labels=allowed_labels, 
                        TREE_FILE=args.TREE_FILE, 
                        output="pruned_tree.nwk")
        if args.ROOT_TREE:
            codes.ROOT_TREE(TREE_FILE="pruned_tree.nwk", output= "rooted_tree.nwk")
