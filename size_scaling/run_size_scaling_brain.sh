snakemake --slurm --default-resources slurm_partition=pinky --set-resources learnMSA:slurm_partition=vision --resources load=10 --jobs 200 --keep-going --rerun-incomplete 
