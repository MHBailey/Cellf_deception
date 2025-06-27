# Cellf_deception
A random forest classifier using gene pairs to assess the tissue of origin

This tool piggybacks off the library and strategy proposed by the multiClassPairs as a way to look at the relationship between features as a classifier in itself. See https://github.com/NourMarzouka/multiclassPairs for the full repo. Here we apply this strategy to many scRNA-seq and bulk RNA-sequencing data of many cell types. 

#Citation 
On bioRxiv soon. 

#Installation (note: I found these libraries do not play well with Apple Silicon M3 chips, especially bioconda::switchBox. It is recommended to run this on a Linux system.) 
1. Install [miniconda](https://www.anaconda.com/docs/getting-started/miniconda/install) 
2. Create an environment `conda create -n cellf`
3. Activate the environment `conda activate cellf`
4. Install the correct libraries.
   conda install r-multiclasspairs 
   conda install r-tidyverse
   conda install bioconda::Biobase
   conda install bioconda::switchBox
   conda install bioconda::snakemake
5. git this repo `git clone https://github.com/MHBailey/Cellf_deception.git`

#Running with SLURM

`pip install snakemake-executor-plugin-slurm`
`snakemake --executor slurm -p all -j 3 --default-resources --slurm-no-account`

#Running without SLURM 

`snakemake -p all -c1` 

#The output of these scripts will provide you with the data cleaning and figures necessary to reproduce our findings. 



