# snakeMSA
A Snakemake pipeline to run easily reproducible benchmarks of various MSA tools.

Runs different MSA tools on the HomFam benchmark. 

## Installation of Snakemake and all MSA tools 
Assuming that conda is already installed, run this once:

`chmod +x ./install.sh ; ./install.sh`

## Usage

Type `conda activate snakeMSA` to activate the conda environment.

Run the pipeline on a workstation:

`snakemake`

Run the pipeline via slurm:

`snakemake --slurm --default-resources slurm_partition=<CPU partition> --set-resources learnMSA:slurm_partition=<GPU partition> --jobs 20 --keep-going`
