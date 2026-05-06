"""
conda create -n DnDs_env python=3
conda activate DnDs_env
conda install bioconda::paml
conda install bioconda::hyphy
conda install conda-forge::biopython
conda install conda-forge::ete3
ete3 build check
cd DnDsCODes
pip install -e .
"""


from Bio import SeqIO
from Bio import AlignIO
from ete3 import Tree
import datetime

def REMOVE_REDUNDANCY(DICT_FASTA: dict =None):
    DICT_TMP = {}
    for ID, seq in DICT_FASTA.items():
        if seq in DICT_TMP:
                print(f"The sequence from ID {ID} is duplicated.")
                print("Just the oldest sequence will be kept.")
                date_previous = DICT_TMP[seq].split("|")[-1]
                date_current = ID.split("|")[-1]
                if date_current < date_previous:
                    DICT_TMP[seq]=ID

        elif len(seq) != 0:
            DICT_TMP[seq]=ID
        else:
            print(f"The ID {ID} was not included in the FASTA dictionary. The sequence has lenght 0.")
        
    DICT_FASTA_NR={}

    for s, k in DICT_TMP.items():
        DICT_FASTA_NR[k] = s

    return DICT_FASTA_NR

def READ_FASTA(FASTA_FILE: str =None, NR: bool = True):
    DICT_TMP={}
    
    for s in SeqIO.parse(FASTA_FILE, "fasta"):
            ID = str(s.description)
            seq=str(s.seq).upper()
            #if "*" not in seq and "x" not in str.lower(seq):
            seq = seq.replace("*","")

            if len(seq) != 0:
                DICT_TMP[ID]=seq

            else:
                print(f"The ID {ID} was not included in the FASTA dictionary. The sequence has lenght 0.")

    if NR:
        return REMOVE_REDUNDANCY(DICT_FASTA=DICT_TMP)
    
    return DICT_TMP
    
def REMOVE_STOP(DICT_FASTA: dict = None):
    STOP_CODONS = ['TAA', 'TAG', 'TGA']
    last_codons = set([ s[-3:].upper() for s in DICT_FASTA.values()])

    thereis = any(elem in last_codons for elem in STOP_CODONS)
    if thereis:
        FINAL_DICT = {k: s[:-3] for k, s in DICT_FASTA.items()}
        return FINAL_DICT
    else:
        return DICT_FASTA

def FastaToPhylip(alignment_file: str = None, output: str = None):
    alignment = AlignIO.read(alignment_file, "fasta")
    with open(output, "w") as out:
        out.write(f"{len(alignment)} {alignment.get_alignment_length()}\n")
        for record in alignment:
            # limita id a 30 chars e garante dois espaços
            name = record.description
            seq  = str(record.seq)
            out.write(f"{name}  {seq}\n")

def PRUNE_TREE(allowed_labels: str = None, TREE_FILE: str = None, output: str = "pruned_tree.nwk"):
    
    tree = Tree(TREE_FILE, format=1)
    tree.prune(allowed_labels, preserve_branch_length=True)

    with open(output, "w") as f:
        f.write(tree.write(format=1))

    return tree

def ROOT_TREE(TREE_FILE: str = None, output: str = "rooted_tree.nwk"):
    tree = Tree(TREE_FILE, format=1)

    def date_extract(tip):
        try:
            date_str = tip.split("|")[-1]
            return datetime.datetime.strptime(date_str, "%Y-%m-%d")
        except ValueError as e:
            print(f"It was not able to extract date from {date_str}. Error {e}")
            return None

    oldest_date = datetime.datetime.max
    oldest_tip = None

    for tip in tree.iter_leaves():
        date = date_extract(tip.name)
        if date and date < oldest_date:
            oldest_date = date
            oldest_tip = tip
    

    if oldest_tip:
        tree.set_outgroup(oldest_tip)
        print(f"Rooted at: {oldest_tip.name}")
    else:
        print("Not found tip to root")
    
    with open(output, "w") as f:
        f.write(tree.write(format=1))
    
    return tree

