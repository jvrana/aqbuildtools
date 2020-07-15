poetry install
poetry lock
poetry build
conda activate aqbt
pip install dist/aqbt*whl
conda env export > environment.yml
