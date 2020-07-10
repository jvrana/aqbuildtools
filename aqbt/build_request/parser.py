import pandas as pd
from typing import *
from itertools import tee
from functools import reduce
import operator
import re

def match_fxn(match: Union[Callable, str, re.Pattern]) -> Callable[[str], bool]:
    if callable(match):
        _match_fxn = match
    else:
        def _match_fxn(x):
            if isinstance(match, re.Pattern):
                return match.match(x)
            else:
                return match == x
    return _match_fxn


def match_value(match: Union[Callable, str, re.Pattern], x: str) -> bool:
    fxn = match_fxn(match)
    return fxn(x)


def find_value(values, match: Union[Callable, str, re.Pattern],
               only_rows: Optional[List[int]] = None,
               only_cols: Optional[List[int]] = None) -> List[Tuple[int, int]]:
    """Return the index"""
    found = []

    for r, row in enumerate(values):
        if not only_rows or r in only_rows:
            for c, val in enumerate(row):
                if not only_cols or c in only_cols:
                    if match_value(match, str(val)):
                        found.append((r, c))
    return found


class Slicer(object):

    def __getitem__(self, x):
        return x


def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)


slicer = Slicer()


def _partition(x, indices: Union[None, int]):
    for a, b in pairwise(indices):
        yield x[a:b]


def _is_empty(row):
    return sum([len(str(r).strip()) for r in row]) == 0


def remove_empty(values):
    return [row for row in values if not _is_empty(row)]


def extract_meta(values):
    r1 = find_value(values, 'Basic DNA Parts')[0][0]
    return remove_empty(values[:r1])


def extract_basic_parts(values):
    r1 = find_value(values, 'Basic DNA Parts')[0][0]
    r2 = find_value(values, 'Composite DNA Parts')[0][0]
    return remove_empty(values[r1:r2])


def extract_composite_parts(values):
    r1 = find_value(values, 'Composite DNA Parts')[0][0]
    return remove_empty(values[r1:])


# parse composite parts


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


def is_square(values):
    """Check that 2D list is square"""
    row_sizes = [len(r) for r in values]
    return len(set(row_sizes)) == 1


def shape(values: List[List]):
    """Return shape of a 2D list"""
    if not is_square(values):
        raise ValueError("Values must be square")
    return len(values), len(values[0])


def empty(s1, s2):
    """Create empty 2D list"""
    new_values = []
    for i in range(s1):
        new_values.append([])
        for j in range(s2):
            new_values[-1].append(None)
    return new_values


def transpose(values):
    num_rows = max([len(r) for r in values])
    transposed = [[] for _ in range(num_rows)]
    for r in range(num_rows):
        for row in values:
            if r < len(row):
                transposed[r].append(row[r])
    return transposed


# def sq_transpose(values: List[List], non_square_ok: bool = False):
#     """Transpose a 2D list"""
#     if not non_square_ok:
#         if not is_square(values):
#             raise ValueError("Values must be square")
#     s = shape(values)
#     new_values = empty(s[1], s[0])
#     for i in range(s[0]):
#         for j in range(s[1]):
#             new_values[j][i] = values[i][j]
#     return new_values

class BuildRequestParsingException(Exception):
    """Generic parsing exception."""


class LocationContext(object):
    """Context manager that relays the location of the parse error."""

    def __init__(self, row, col):
        self.row = row
        self.col = col

    def __enter__(self):
        pass

    def __exit__(self, exc_type, exc_value, traceback):
        if exc_value:
            raise BuildRequestParsingException(
                "{} near ({}, {}).\n{}".format(exc_type, self.row, self.col, exc_value))


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


# is_square([[]])

values = to_cell_values(data['values'])
parse_composite_parts(values);

part_json = pd.DataFrame(extract_basic_parts(values)[2:],
                         columns=extract_basic_parts(values)[1]).T.to_dict()
part_json

# validate