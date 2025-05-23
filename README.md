# snakeMSA

A Snakemake pipeline to run easily reproducible benchmarks of various MSA tools. It can be used to reproduce the results of our learnMSA2 paper [1].

The pipeline runs different aligners on a data collection and scores the results using reference alignments. Each aligner is installed and run in an isolated environment. The aligner parameters are defined in a unified configuration file.


## Installation of Snakemake and all MSA tools 

Assuming that conda is already installed, run this once:

`chmod +x ./install.sh ; ./install.sh`


## Aquire the HomFam collection

Run

`chmod +x ./data/homfam/get_homfam.sh ; ./data/homfam/get_homfam.sh`

to download and preprocess the HomFam collection.


## Run on custom data

The pipeline can be run on any data colletion. In your configuration file (for example see `config/default.json`) define the field `data_path` to point to a directory with 2 folders: `aligned` and `unaligned`. Both directories should contain fasta files with matching filenames. The aligned sequences can be a subset of the unaligned sequences, but the sequence identifiers must match. The unaligned sequences will be aligned by the tools specified in the configuration file and the aligned references will be used to score the output.


## Run locally

To align HomFam with a variety of aligners run

`conda run -n snakeMSA_base_env snakemake --cores all --configfile configs/benchmark_all.json`.

To align HomFam with the 3 variants of learnMSA used in the learnMSA2 paper run

`conda run -n snakeMSA_base_env snakemake --cores all --configfile configs/benchmark_learnMSA.json`.


## Run on a SLURM cluster

`snakemake`

Run the pipeline via slurm:

`conda run -n snakeMSA_base_env snakemake --executor slurm --default-resources slurm_partition=<your partition> --slurm-logdir slurm --jobs 20 --configfile configs/benchmark_learnMSA.json`.

By default, learnMSA expects a GPU partition. To run learnMSA on a CPU partition (not recommended), change the `use_gpu` field in `benchmark_learnMSA.json`.



[1] Becker, F., & Stanke, M. (2024). learnMSA2: deep protein multiple alignments with large language and hidden Markov models. Bioinformatics, 40(Supplement_2), ii79-ii86.
