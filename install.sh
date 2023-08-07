conda install -n base -c conda-forge mamba
mamba env create --name snakeMSA --file environment.yaml
mamba activate snakeMSA
pip install magus-msa