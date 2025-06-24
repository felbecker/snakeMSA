# snakeMSA

A Snakemake pipeline to run easily reproducible benchmarks of various MSA tools. It can be used to reproduce the results of our learnMSA2 paper [1].

The pipeline runs different aligners on a data collection and scores the results using reference alignments. 
Each aligner is installed and run in an isolated environment. The aligner parameters are defined in a unified configuration file.


## Installation of Snakemake and all MSA tools 

Assuming that conda is already installed, run this once:

`chmod +x ./install.sh ; ./install.sh`


## Aquire the HomFam collection

The HomFam collection can be downloaded via:

`chmod +x ./data/homfam/get_homfam.sh ; ./data/homfam/get_homfam.sh`

The extended HomFam collection [2] can be downloaded via:

`chmod +x ./data/ext_homfam/get_ext_homfam.sh ; ./data/ext_homfam/get_ext_homfam.sh`


## Run on custom data

The pipeline can be run on any data colletion. In your configuration file (for example see `config/default.json`) define the field `data_path` to point to a directory with 2 folders: `aligned` and `unaligned`. 
Both directories should contain fasta files with matching filenames (and extensions). 
The aligned sequences can be a subset of the unaligned sequences, but the sequence identifiers must match. 
The unaligned sequences will be aligned by the tools specified in the configuration file and the aligned sequences will be used as references to score the outputs.


## Run locally

To align HomFam with a variety of aligners run:

`conda run -n snakeMSA_base_env snakemake --cores all --configfile configs/all.json`.

To align HomFam with the 3 variants of learnMSA used in the learnMSA2 paper run:

`conda run -n snakeMSA_base_env snakemake --cores all --configfile configs/learnMSA.json`.

There are also `*_huge.json` variants available to run the large scale experiments.


## Run on a SLURM cluster

The [Snakemake SLURM executor plugin](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html) can be used to run alignment jobs on a SLURM cluster.

The recommended way is to provide a profile named `config.yaml` with the following minimal content:

```
executor: slurm
jobs: 30
max-jobs-per-second: 10
default-resources:
 - slurm_partition=<your CPU or GPU partition>
slurm-logdir: slurm
```

Run the pipeline via slurm:

`conda run -n snakeMSA_base_env snakemake --profile <path-to-folder-containing-config.yaml> --configfile configs/learnMSA.json`.

Please note: by default, learnMSA expects a GPU partition (and we recommend to use one!). For CPU change `use_gpu` field in `learnMSA.json` to False.


## Customizing job resources

Job resources can be customized in the configuration file. See `configs/default.json`.


[1] Becker, F., & Stanke, M. (2024). learnMSA2: deep protein multiple alignments with large language and hidden Markov models. Bioinformatics, 40(Supplement_2), ii79-ii86.
[2] Deorowicz, S., Debudaj-Grabysz, A., & Gudyś, A. (2016). FAMSA: Fast and accurate multiple sequence alignment of huge protein families. Scientific reports, 6(1), 33964.
