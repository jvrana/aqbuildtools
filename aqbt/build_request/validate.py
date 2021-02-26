import json
from os.path import abspath
from os.path import dirname
from os.path import join
from typing import Dict
from typing import List
from typing import Union

import jsonschema

from .exceptions import parsing_location

here = dirname(abspath(__file__))

Part = Dict

with open(join(here, "schema.json")) as f:
    schema = json.load(f)


class BuildRequestParsingError(Exception):
    """Generic build request parsing error."""


def schema_validate(x: Dict[str, Union[Dict, str, int, float]]) -> None:
    return jsonschema.validate(x, schema)


def validate_part(part: Part):
    try:
        schema_validate(part)
    except jsonschema.exceptions.ValidationError as e:
        raise BuildRequestParsingError(str(e)) from e

    if part["partType"] == "basic part":
        if "length" in part and not part["length"] == len(part["sequence"]):
            raise BuildRequestParsingError(
                "Part {}: 'length' ({}) does not match the length"
                " of the provided sequence ({})".format(
                    part["name"], part["length"], len(part["sequence"])
                )
            )

    if part["name"] == "nan":
        with parsing_location(part["name"]):
            raise BuildRequestParsingError(
                '"nan" is an invalid part name for a basic part'
            )

    if part["partType"] == "composite part":
        for _part in part["parts"]:
            with parsing_location(_part):
                if _part == "nan":
                    raise BuildRequestParsingError(
                        '"nan" is an invalid part name for a sub part'
                    )


def validate_part_list(parts: List[Part], fast_fail: bool = True):
    for part in parts:
        with parsing_location(part["name"]):
            validate_part(part)

    composite_parts = [part for part in parts if part["partType"] == "composite part"]
    basic_parts = [part for part in parts if part["partType"] == "basic part"]

    # {name: part} dictionary
    part_dict = {}

    # check for name conflicts
    for part in basic_parts + composite_parts:
        if part["name"].lower() in part_dict:
            with parsing_location(part["name"]):
                raise BuildRequestParsingError(
                    "Part name conflict for {}".format(part["name"])
                )
        else:
            part_dict[part["name"].lower()] = part

    # check for references
    for composite_part in composite_parts:
        for sub_part_name in composite_part["parts"]:
            if sub_part_name.lower() not in part_dict:
                with parsing_location(composite_part["name"]):
                    raise BuildRequestParsingError(
                        "Subpart '{}.{}' is missing a definition".format(
                            composite_part["name"], sub_part_name
                        )
                    )
