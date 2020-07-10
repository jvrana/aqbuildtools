"""
Methods for parsing a build request document.

The primary input for this is a List[List[str]] representing the
excel or csv sheet. This can be obtained by parsing an excel sheet, using
the google sheet API, etc.
"""

import pandas as pd
from typing import TypeVar, Union, Callable, List, Iterable, Tuple, Optional
from itertools import tee
import re
from .cell_value import CellValue
from .exceptions import LocationContext, BuildRequestParsingException

T = TypeVar('T')

##########################
# Utilities
##########################

def _dispatch_match_fxn(match: Union[Callable, str, re.Pattern]) -> Callable[[str], bool]:
    """Generate a match function from a string, Callable, or regex pattern."""
    if callable(match):
        _match_fxn = match
    else:
        def _match_fxn(x):
            if isinstance(match, re.Pattern):
                return match.match(x)
            else:
                return match == x
    return _match_fxn


def match_value(match: Union[Callable, str, re.Pattern], value: str) -> bool:
    """Match a string, Callable, or regex pattern with a value"""
    fxn = _dispatch_match_fxn(match)
    return fxn(value)


def find_value(values, match: Union[Callable, str, re.Pattern],
               only_rows: Optional[List[int]] = None,
               only_cols: Optional[List[int]] = None) -> List[Tuple[int, int]]:
    """Return the indices [(int, int), ...] of the cells that match the provided
    matching string, function, regex pattern."""
    found = []

    for r, row in enumerate(values):
        if not only_rows or r in only_rows:
            for c, val in enumerate(row):
                if not only_cols or c in only_cols:
                    if match_value(match, str(val)):
                        found.append((r, c))
    return found

def pairwise(iterable: Iterable[T]) -> Iterable[Tuple[T, T]]:
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)


def _partition(x: List, indices: Union[None, int]) -> List:
    """Partition a list using a list of indices."""
    for a, b in pairwise(indices):
        yield x[a:b]


def _is_empty(row: List[str]) -> bool:
    """check that a row of string values is empty"""
    return sum([len(str(r).strip()) for r in row]) == 0


def remove_empty(values: List[List[str]]) -> List[List[str]]:
    """Remove empty rows"""
    return [row for row in values if not _is_empty(row)]

def is_square(values):
    """Check that 2D list is square"""
    row_sizes = [len(r) for r in values]
    return len(set(row_sizes)) == 1

def transpose(values):
    num_rows = max([len(r) for r in values])
    transposed = [[] for _ in range(num_rows)]
    for r in range(num_rows):
        for row in values:
            if r < len(row):
                transposed[r].append(row[r])
    return transposed

##########################
# Extract sections of
# the build request values
##########################

def extract_meta(values: List[List[CellValue]]) -> List[List[CellValue]]:
    r1 = find_value(values, 'Basic DNA Parts')[0][0]
    return remove_empty(values[:r1])


def extract_basic_parts(values):
    r1 = find_value(values, 'Basic DNA Parts')[0][0]
    r2 = find_value(values, 'Composite DNA Parts')[0][0]
    return remove_empty(values[r1:r2])


def extract_composite_parts(values):
    r1 = find_value(values, 'Composite DNA Parts')[0][0]
    return remove_empty(values[r1:])


########################
# Parse extracted values
########################


def values_to_json(values, key_column: int, fill_empty_key_from_above: bool = False):
    data = {}
    for row in values:
        key = row[key_column].strip()
        if fill_empty_key_from_above and key == '':
            key = prev_key
        data.setdefault(key, list())
        data[key].append(row[key_column + 1:])
        prev_key = key
    return data

def parse_composite_parts(values):
    composite = extract_composite_parts(values)

    # locate indices of new "Collections"
    indices = find_value(composite, re.compile('Collection Name'))
    rows = [None] + [i[0] for i in indices] + [None]

    # partition the values on the indices
    partitioned = list((_partition(composite, rows)))[1:]
    parsed_json_arr = []
    for p in partitioned:
        with LocationContext(*p[0][0].rc()):
            parsed_json_arr.append(values_to_json(p, 0, True))

    def _make_part_list(parsed_json):
        with LocationContext(*parsed_json['Name:'][0][0].rc()):
            names = transpose(parsed_json['Name:'])
        with LocationContext(*parsed_json['Parts:'][0][0].rc()):
            parts = transpose(parsed_json['Parts:'])
        with LocationContext(*parsed_json['Description:'][0][0].rc()):
            descriptions = transpose(parsed_json['Description:'])
        collection_name = parsed_json['Collection Name:'][0][0]

        part_list = []
        for _name, _description, _parts in list(zip(names, descriptions, parts)):
            part_list.append({
                'name': _name[0],
                'collection': collection_name,
                'type': 'composite part',
                'description': _description[0],
                'parts': _parts
            })
        return part_list

    part_list = []
    for p in parsed_json_arr:
        part_list += _make_part_list(p)
    return part_list


values = to_cell_values(data['values'])
parse_composite_parts(values);

part_json = pd.DataFrame(extract_basic_parts(values)[2:],
                         columns=extract_basic_parts(values)[1]).T.to_dict()
part_json

# validate