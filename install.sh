#!/bin/bash

# iterate over the files *_env.yaml in envs
# for each strip .yaml and create a conda environment with the name
# and the yaml file as requirements
for file in envs/*.yaml; do
    # strip the .yaml and the path
    name=$(basename "$file" .yaml)
    echo "Creating conda environment $name from $file"
    # Check if the environment already exists
    if conda info --envs | grep -q "^$name\s"; then
        echo "Environment $name already exists. Skipping."
    else
        conda env create --name "$name" --file "$file" || {
            echo "Failed to create environment $name. Exiting."
            exit 1
        }
    fi
done

# additional install required to align with pLM support
conda activate learnMSA_env && pip install torch==2.2.1 tf-keras==2.17.0
conda deactivate