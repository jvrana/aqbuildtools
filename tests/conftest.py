from os.path import abspath
from os.path import dirname
from os.path import join

import pytest
import toml
from pydent import AqSession

from aqbt import AquariumBuildTools
from aqbt.aquarium.registry import LabDNARegistry
from aqbt.tools import config_to_sessions
from aqbt.tools import parse_config
from benchlingapi import Session as BenchlingSession


here = abspath(dirname(__file__))


@pytest.fixture(scope="session")
def config_path():
    return join(here, "secrets", "test_config.toml")


@pytest.fixture(scope="session")
def config_alt_path():
    return join(here, "secrets", "test_config.alt.toml")


@pytest.fixture(scope="session")
def config(config_path):
    with open(config_path) as f:
        return parse_config(toml.load(f))


@pytest.fixture(scope="session")
def sessions(config):
    return config_to_sessions(config)


@pytest.fixture(scope="session")
def fixtures_path():
    return join(here, "fixtures")


@pytest.fixture(scope="session")
def registry(sessions) -> LabDNARegistry:
    return sessions["default"]["registry"]


@pytest.fixture(scope="session")
def benchling(sessions) -> BenchlingSession:
    return sessions["default"]["benchling"]


@pytest.fixture(scope="session")
def aquarium(sessions) -> AqSession:
    return sessions["default"]["aquarium"]


@pytest.fixture(scope="session")
def tools(config):
    return AquariumBuildTools(config)
