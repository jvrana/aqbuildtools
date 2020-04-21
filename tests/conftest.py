import pytest
from os.path import abspath, dirname, join
from aqbt.cli import parse_config, config_to_sessions
import toml

here = abspath(dirname(__file__))


@pytest.fixture(scope='session')
def config_path():
    return join(here, 'secrets', 'test_config.toml')


@pytest.fixture(scope='session')
def sessions(config_path):
    with open(config_path, 'r') as f:
        return config_to_sessions(parse_config(toml.load(f)))

@pytest.fixture(scope='session')
def fixtures_path():
    return join(here, 'fixtures')