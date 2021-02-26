import json
from os.path import abspath
from os.path import dirname
from os.path import join
from typing import Dict
from typing import Union

import jsonschema

here = dirname(abspath(__file__))

with open(join(here, "schema.json")) as f:
    schema = json.load(f)


def schema_validate(x: Dict[str, Union[Dict, str, int, float]]) -> None:
    return jsonschema.validate(x, schema)
