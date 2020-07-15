# updates the conda environment for the aqbt
conda activate aqbt
pip install -r requirements.txt
conda env export > conda_environment.yml
conda env update --prefix ./env --file conda_environment.yml --prune
conda deactivate aqbt