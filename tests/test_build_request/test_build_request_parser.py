import pytest
import pandas as pd
from os.path import join
from aqbt.build_request.parser import parse_composite_parts, parse_basic_parts, parse_parts, CellValue
from aqbt.build_request.validate import validate_part, validate_part_list
from aqbt.build_request.exceptions import LocationContext


@pytest.fixture()
def values(fixtures_path):
    df = pd.read_csv(join(fixtures_path, 'build_request_fixtures', 'yg-designs.csv'))
    values = df.values.tolist()
    values = CellValue.to_cell_values(values)
    return values

def test_parse_and_validate_composite_parts(values):
    parsed = parse_composite_parts(values)
    assert parsed
    for p in parsed:
        with LocationContext(*p['name'].rc()):
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