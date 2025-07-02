configfile: 'config.yaml'

import datetime as dt
from pathlib import Path
import os


DATE = dt.datetime.today().strftime("%Y_%m_%d")
DAY = dt.datetime.today().strftime("%d")

localrules: all

rule agg_data:
    input:
        code=config["Clean"]
    output:
        outlast="Processed_data/cleandata_combined_training_line.MHB.csv",
        outcd ="Processed_data/cleandata_combined_training_tissue.MHB.csv",
    resources:
        runtime='210m',
        mem_mb=30000
    shell:
        '''
        mkdir -p Figures Processed_data
        Rscript --vanilla --quiet {input.code}
        '''

rule primary_cell_classify:
    input:
        code=config["Primary"],
        newin="Processed_data/cleandata_combined_training_tissue.MHB.csv"
    output:
        out="Figures/Primary_chandra_plot.pdf"
    resources:
        runtime='210m',
        mem_mb=30000
    shell:
        '''
        Rscript --vanilla --quiet {input.code}
        '''

rule depmap:
    input:
        code=config["Depmap"],
        newin="Processed_data/cleandata_combined_training_line.MHB.csv"
    output:
        out="Figures/Depmap_hmc3_plot.pdf"
    resources:
        runtime='210m',
        mem_mb=30000
    shell:
       '''
       Rscript --vanilla --quiet {input.code} 
       '''

rule all:
    input:
        figPrimary="Figures/Primary_chandra_plot.pdf",
        figDM="Figures/Depmap_hmc3_plot.pdf"
       
