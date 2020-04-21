init:
	pip install pip -U
	curl -sSL https://raw.githubusercontent.com/sdispater/poetry/master/get-poetry.py | python
	poetry install
	poetry run pre-commit install


clean:
	rm -rf dist
	rm -rf pip-wheel-metadata
	rm -rf docs
	rm -rf .pytest_cache


test:
	poetry run python -m pytest


lint:
	poetry run pylint -E pydent


docs:
	echo "No documentation"


format:
	poetry run keats version up
	poetry run black lobio tests


lock:
	poetry run keats version up
	poetry update


build:
	poetry run keats version up
	poetry build


release:
	poetry run keats run release


klocs:
	find . -name '*.py' | xargs wc -l
