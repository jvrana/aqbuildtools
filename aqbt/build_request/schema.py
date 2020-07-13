from os.path import join, abspath, dirname
import json
import jsonschema
from typing import Union, Dict

here = dirname(abspath(__file__))

with open(join(here, "schema.json"), "r") as f:
    schema = json.load(f)


def schema_validate(x: Dict[str, Union[Dict, str, int, float]]) -> None:
    return jsonschema.validate(x, schema)
