from typing import Dict
from typing import Union

import toml
from benchlingapi import Session as BenchlingSession
from pydent import AqSession

from aqbt.aquarium.registry import LabDNARegistry
from aqbt.aquarium.registry import BenchlingRegistryConnector


class ConfigParsingException(Exception):
    pass


def parse_config(config: Dict[str, Union[str, dict, int]]) -> Dict[str, Dict[str, str]]:
    sessions = {}

    # validation checks
    errors = []
    if "session" not in config:
        errors.append("Config must have at least one '[session.<NAME>]' key.")

    if "aquarium" not in config:
        errors.append("Config must have at least one '[aquarium.<SCOPE>] key.")

    if "benchling" not in config:
        errors.append("Config must have at least one '[benchling.<SCOPE>] key.")

    if not errors:
        for session_name, session_info in config["session"].items():

            aquarium_scope = session_info["aquarium"]
            benchling_scope = session_info["benchling"]

            if aquarium_scope not in config["aquarium"]:
                errors.append(
                    "Aquarium is missing a definition for the [aquarium.{scope}] "
                    "scope found in the [session.{session}] session definition.".format(
                        scope=aquarium_scope, session=session_name
                    )
                )

            if benchling_scope not in config["benchling"]:
                errors.append(
                    "Benchling is missing a definition for the [aquarium.{scope}] "
                    "scope found in the [session.{session}] session definition.".format(
                        scope=benchling_scope, session=session_name
                    )
                )

    if errors:
        raise ConfigParsingException(
            "Could not parse config due to the following errors:\n{}".format(
                "\n".join(
                    "({i}) - {err}".format(i=i, err=err) for i, err in enumerate(errors)
                )
            )
        )

    for session_name, session_info in config["session"].items():
        sessions[session_name] = {
            "aquarium": config["aquarium"][aquarium_scope],
            "benchling": config["benchling"][benchling_scope],
        }
    return sessions


def config_to_sessions(config: Dict[str, Dict[str, str]]):
    sessions = {}
    for session_name, session_info in config.items():
        sessions[session_name] = {"aquarium": None, "benchling": None, "registry": None}

        sessions[session_name]["aquarium"] = AqSession(
            login=session_info["aquarium"]["login"],
            password=session_info["aquarium"]["password"],
            aquarium_url=session_info["aquarium"]["url"],
        )

        sessions[session_name]["benchling"] = BenchlingSession(
            api_key=session_info["benchling"]["apikey"]
        )

        registry_connector = BenchlingRegistryConnector(
            api=sessions[session_name]["benchling"],
            initials=session_info["benchling"]["initials"],
            schema=session_info["benchling"]["schema"],
            prefix=session_info["benchling"]["prefix"],
            folder_id=session_info["benchling"]["folder_id"],
            registry_id=session_info["benchling"]["id"],
        )  #: klavins lab Benchling registry connector

        sessions[session_name]["registry"] = LabDNARegistry(
            connector=registry_connector, aqsession=sessions[session_name]["aquarium"],
        )  #: klavins lab Benchling registry
    return sessions


def generate_example_config():
    return {
        "aquarium": {"default": {"login": "", "password": "", "url": ""}},
        "benchling": {
            "default": {
                "apikey": "",
                "initials": "",
                "schema": "",
                "prefix": "",
                "folder_id": "",
                "id": "",
            }
        },
        "session": {"default": {"aquarium": "default", "benchling": "default"}},
    }


from aqbt.aquarium.genome_builder import aq_to_gff


class AquariumBuildTools:
    def __init__(self, config: Dict[str, Dict[str, str]]):
        self.sessions = config_to_sessions(config)

    @classmethod
    def from_toml(cls, path: str):
        with open(path, "r") as f:
            config_dict = toml.load(f)
        config = parse_config(config_dict)
        return cls(config)
