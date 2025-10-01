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

`cd data/ext_homfam && chmod +x ./get_ext_homfam.sh && ./get_ext_homfam.sh`


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

## Specific experiments

### Size scaling

This experiment explores how alignment accuracy behaves relative to sequence count.
We provide a script called `make_nested_sets.py` that takes an existing benchmark collection as well as a list of 
sequence counts `L`.
It will only consider unaligned datasets in the benchmark collection with have at least as many sequences as the maximum number in the provided list.
It also assumes that the number of sequences in each reference alignment is at most the minimum number in the provided list.
The script will randomly sample `n` series of nested datasets per dataset in the original collection with the following properties:
- Dataset `i` in the series has `L[i]` sequences.
- All datasets in a series contain the reference sequences.
- Each datasets in a series contains the preceeding dataset with less sequences as a subset.

To generate data for the size scaling experiment from the learnMSA2 paper, run:

`conda run -n learnMSA_env python tools/make_nested_sets.py --benchmark data/ext_homfam/ext_homfam_huge --length 1e3 1e4 5e4 1e5 2.5e5 --n_series 4 --prefix data/size_scaling_homfam_huge`

Use the regular pipeline to align `data/size_scaling_homfam_huge`.

Then run the pipeline on the generated collection.

## References

[1] Becker, F., & Stanke, M. (2024). learnMSA2: deep protein multiple alignments with large language and hidden Markov models. Bioinformatics, 40(Supplement_2), ii79-ii86.

[2] Deorowicz, S., Debudaj-Grabysz, A., & Gudyś, A. (2016). FAMSA: Fast and accurate multiple sequence alignment of huge protein families. Scientific reports, 6(1), 33964.
