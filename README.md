conda create -n DnDs_env python=3  
conda activate DnDs_env
conda install bioconda::paml
conda install bioconda::hyphy
conda install conda-forge::biopython
conda install conda-forge::ete3
ete3 build check
cd DnDsCODes
pip install -e .
