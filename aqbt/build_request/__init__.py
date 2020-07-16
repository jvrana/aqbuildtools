"""
BuildRequest document parser
"""

from .parser import parse_composite_parts, parse_basic_parts, parse_parts, CellValue
from .schema import schema, schema_validate
from os.path import exists

import pandas as pd


def parse(data, format):
    is_file = False
    if isinstance(data, str):
        if exists(data):
            is_file = True
    if format == 'csv':
        pd.read_csv()

