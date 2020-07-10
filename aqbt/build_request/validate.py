from os.path import join, abspath, dirname
import json
import jsonschema
from typing import Union, Dict, List
here = dirname(abspath(__file__))

Part = Dict

with open(join(here, 'schema.json'), 'r') as f:
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

    if part['partType'] == 'basic part':
        if 'length' in part and not part['length'] == len(part['sequence']):
            raise BuildRequestParsingError("Part {}: 'length' ({}) does not match the length"
                                           " of the provided sequence ({})".format(
                part['name'], part['length'], len(part['sequence'])
            ))


def validate_part_list(parts: List[Part]):
    for part in parts:
        validate_part(part)

    composite_parts = [part for part in parts if part['partType'] == 'composite part']
    basic_parts = [part for part in parts if part['partType'] == 'basic part']

    part_dict = {}

    # check for name conflicts
    for part in basic_parts + composite_parts:
        if part['name'] in part_dict:
            raise BuildRequestParsingError("Part name conflict for {}".format(part['name']))

    # check for references
    for composite_part in composite_parts:
        for sub_part_name in composite_part['parts']:
            if sub_part_name not in part_dict:
                raise BuildRequestParsingError("Subpart '{}.{}' is missing a definition".format(
                    composite_part['name'], sub_part_name
                ))
