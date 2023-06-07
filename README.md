# snakeMSA
A Snakemake pipeline to run easily reproducible benchmarks of various MSA tools.

Runs different MSA tools on the extHomFam 2 benchmark which consists of 5 sub-folders (sorted by MSA depth). Outputs {tool}.tbl files with SP, modeler, TC and column scores per protein family. 

## Installation of Snakemake and all MSA tools 
Assuming that conda is already installed, run this once:

`chmod +x ./install.sh ; ./install.sh`

## Usage

Type `conda activate snakeMSA` to activate the conda environment.

Run the pipeline on a workstation:

`snakemake`

Run the pipeline via slurm:

`snakemake --slurm --default-resources slurm_partition=<CPU partition> --set-resources learnMSA:slurm_partition=<GPU partition> --jobs 50`
