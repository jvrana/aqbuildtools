# updates the conda environment for the aqbt
conda activate aqbt
pip install $(ls dist/aqbt*whl)
conda env export > conda_environment.yml
conda env update --prefix ./env --file conda_environment.yml --prune
conda deactivate aqbt