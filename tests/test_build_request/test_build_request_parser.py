import pytest
import pandas as pd
from os.path import join
from aqbt.build_request.parser import parse_composite_parts, CellValue


def test_parse_composite_parts(fixtures_path):
    df = pd.read_csv(join(fixtures_path, 'yg-designs.csv'))
    values = df.values.tolist()
    values = CellValue.to_cell_values(values)
    parsed = parse_composite_parts(values)
    for p in parsed[:10]:
        print(p)
