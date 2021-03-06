from os.path import join

import pandas as pd
import pytest

from aqbt.build_request.exceptions import LocationContext
from aqbt.build_request.parser import CellValue
from aqbt.build_request.parser import parse_basic_parts
from aqbt.build_request.parser import parse_composite_parts
from aqbt.build_request.parser import parse_parts
from aqbt.build_request.parser import row_is_empty
from aqbt.build_request.validate import validate_part
from aqbt.build_request.validate import validate_part_list


@pytest.fixture()
def values(fixtures_path):
    df = pd.read_csv(join(fixtures_path, "build_request_fixtures", "yg-designs.csv"))
    values = df.values.tolist()
    values = CellValue.to_cell_values(values)
    return values


class TestParseUtilities:
    def test_is_empty1(self):
        assert row_is_empty(["nan"] * 10)

    def test_is_empty2(self):
        assert row_is_empty(["nan "] * 10 + [""] + [" "])

    def test_is_not_empty(self):
        assert not row_is_empty(["nan "] * 10 + [""] + [" "] + ["x"])


def test_parse_and_validate_composite_parts(values):
    parsed = parse_composite_parts(values)
    assert parsed
    for p in parsed:
        with LocationContext(*p["name"].rc()):
            validate_part(p)


def test_parse_and_validate_basic_parts(values):
    parsed = parse_basic_parts(values)
    for p in parsed:
        print(p)
        try:
            validate_part(p)
        except Exception as e:
            raise e


def test_parsing(values):
    parts = parse_parts(values)
    validate_part_list(parts)
