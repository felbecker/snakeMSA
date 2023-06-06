# snakeMSA
A Snakemake pipeline to run easily reproducible benchmarks of various MSA tools.

Runs different MSA tools on the extHomFam dataset which consists of 5 sub-folders (sorted by MSA depth) with fasta files. Outputs {tool}.out files with SP, modeler, TC and column scores per protein family in the benchmark. 

## Installation (once, assuming that conda is already installed)
`chmod +x ./install.sh ; ./install.sh`

## Usage

Type `conda activate snakeMSA` to activate the conda environment.

Run the pipeline on a workstation:

`snakemake`

Run the pipeline via slurm:

`snakemake --slurm --default-resources slurm_account=<your SLURM account> slurm_partition=<your SLURM partition>`

Limitations:
Currently CPU-only, so our own tool, [learnMSA](https://github.com/Gaius-Augustus/learnMSA), which uses GPU, has to be run separately. 
Todo: Allow to run learnMSA on a different slurm partition than the rest of the tools.
