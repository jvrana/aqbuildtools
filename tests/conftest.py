from os.path import abspath
from os.path import dirname
from os.path import join

import pytest
import toml

from aqbt.tools import config_to_sessions
from aqbt.tools import parse_config

here = abspath(dirname(__file__))


@pytest.fixture(scope="session")
def config_path():
    return join(here, "secrets", "test_config.toml")


@pytest.fixture(scope="session")
def config(config_path):
    with open(config_path, "r") as f:
        return parse_config(toml.load(f))


@pytest.fixture(scope="session")
def sessions(config):
    return config_to_sessions(config)


@pytest.fixture(scope="session")
def fixtures_path():
    return join(here, "fixtures")


@pytest.fixture(scope="session")
def registry(sessions):
    return sessions["default"]["registry"]


@pytest.fixture(scope="session")
def benchling(sessions):
    return sessions["default"]["benchling"]


@pytest.fixture(scope="session")
def aquarium(sessions):
    return sessions["default"]["aquarium"]
