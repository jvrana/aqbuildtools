# Creates a new conda environment with the necessary aqbt requirments
poetry install
poetry export -f requirements.txt > requirements.txt
conda env create aqbt python=3.8
conda activate aqbt
pip install -r requirements.txt
conda env export > conda_environment.yml
conda deactivate aqbt