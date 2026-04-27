# DnDsCODes 🧬

**DnDsCODes** is a Python package developed to automate and simplify the preprocessing of sequence alignments and phylogenetic trees. It prepares your structural data for natural selection (dN/dS) and genomic epidemiology downstream analyses in software such as PAML and HyPhy.

## 🚀 Features

* **Redundancy Removal:** Filters identical sequences. In case of identical sequences, it prioritizes keeping the oldest sample (based on the collection date extracted from the header).
* **Stop Codon Cleaning:** Detects and removes stop codons at the end of the alignment, preventing *frameshift* errors in dN/dS analyses.
* **Format Conversion:** Easily converts alignments from FASTA to PHYLIP (`.phy`) format.
* **Tree Pruning:** Automatically synchronizes the phylogenetic tree with the cleaned alignment, removing missing leaves/taxa.
* **Temporal Rooting:** Roots the phylogenetic tree using the sequence with the oldest collection date.
* **Split by Year:** Automatically splits large datasets (alignments and trees) into annual subsets for stratified temporal analyses.

---

## ⚙️ Prerequisites and Installation

We highly recommend using `conda` (Miniconda or Anaconda) to manage your environment and dependencies, ensuring no package conflicts occur.

### 1. Creating the virtual environment and installing dependencies
Open your terminal and run the following commands:

```bash
# Create an environment named DnDs_env with Python 3
conda create -n DnDs_env python=3

# Activate the environment
conda activate DnDs_env

# Install evolutionary biology tools (Optional, but recommended for downstream analysis)
conda install bioconda::paml
conda install bioconda::hyphy

# Install DnDsCODes direct dependencies
conda install conda-forge::biopython
conda install conda-forge::ete3

# Verify ete3 installation
ete3 build check
```

### 2. Installing DnDsCODes
After cloning this repository, navigate to the package's main directory and install it in editable mode:

```bash
git clone [https://github.com/YOUR_USERNAME/DnDsCODes.git](https://github.com/YOUR_USERNAME/DnDsCODes.git)
cd DnDsCODes
pip install -e .
```
*(Remember to replace `YOUR_USERNAME` with your actual GitHub username).*

---

## ⚠️ Input Data Format (Important!)

For the temporal rooting (`--ROOT_TREE`) and annual splitting (`--SPLIT_BY_YEAR`) functions to work correctly, **the headers of your FASTA file and the leaves of your Phylogenetic Tree must follow a specific pattern**. 

The collection date (in `YYYY-MM-DD` or `YYYY` format) must be the **last piece of information**, separated by a pipe (`|`).

**Correct example in the `.fasta` file:**
```text
>Isolate_Alpha|Brazil|EpiCoV|2020-05-12
ATGCGTACGTTAG...
>Isolate_Beta|Argentina|EpiCoV|2021-11-03
ATGCGTACGTTAG...
```

---

## 💻 How to Use

Once the environment is activated, you can run the main pipeline using the command line.

### Available Parameters:

| Argument | Required | Description |
| :--- | :---: | :--- |
| `--ALN_FASTA_FILE` | **Yes** | Path to the alignment file in FASTA format. |
| `--TREE_FILE` | **Yes** | Path to the phylogenetic tree file in Newick (`.nwk`) format. |
| `--ROOT_TREE` | No | *Flag*. If passed, roots the tree using the oldest sample. |
| `--SPLIT_BY_YEAR` | No | *Flag*. If passed, splits the alignments and trees by collection year. |
| `--SAVE_PHYLIP` | No | *Flag*. If passed, also saves the resulting alignment in PHYLIP format. |

### Example 1: Basic Usage (Cleaning and Pruning)
Performs redundancy removal, stop codon exclusion, and prunes the phylogenetic tree to match the alignment.

```bash
python -m DnDsCODes --ALN_FASTA_FILE my_alignment.fasta --TREE_FILE my_tree.nwk
```
*Outputs:* * `aln_pruned_withoutStop.fasta` 
* `pruned_tree.nwk`

### Example 2: Full Pipeline
Executes basic cleaning, saves outputs in PHYLIP (for PAML/HyPhy), temporally roots the tree, and splits everything by year.

```bash
python -m DnDsCODes \
  --ALN_FASTA_FILE my_alignment.fasta \
  --TREE_FILE my_tree.nwk \
  --SAVE_PHYLIP \
  --ROOT_TREE \
  --SPLIT_BY_YEAR
```
*Outputs (for each detected year, e.g., 2020, 2021):* * `2020_aln_pruned_withoutStop.fasta`
* `2020_aln_pruned_withoutStop.phy`
* `2020_pruned_tree.nwk`
* `2020_rooted_tree.nwk`

---

## 🤝 Contributing
Contributions are welcome! Feel free to open *Issues* reporting bugs, suggesting new features, or submitting *Pull Requests*.
