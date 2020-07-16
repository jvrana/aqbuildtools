"""
Methods for parsing a build request document.

The primary input for this is a List[List[str]] representing the
excel or csv sheet. This can be obtained by parsing an excel sheet, using
the google sheet API, etc.
"""

import pandas as pd
from typing import TypeVar, Union, Callable, List, Iterable, Tuple, Optional, Dict
from itertools import tee
import re
from .cell_value import CellValue
from .exceptions import parsing_location, BuildRequestParsingException
from .validate import validate_part_list

T = TypeVar("T")


# TODO: parse these keys
class Keys:
    # keys used in the build request document

    Basic_DNA_Parts = "Basic DNA Parts"
    Composite_DNA_Parts = "Composite DNA Parts"
    Collection_Name = "Collection Name"
    Part_Name = "Part Name"
    Description = "Description (Optional)"
    Sequence = "Sequence"
    Length = "length (bp)"
    Source = "Source (Optional)"
    Role = "Role"


##########################
# Utilities
##########################

def _is_nan(x: Union[str, CellValue]) -> bool:
    if x is None:
        return True
    elif x.strip().lower() == 'nan':
        return True
    elif x.strip().lower() == 'None':
        return True
    elif x.strip() == '':
        return True
    return False

def _dispatch_match_fxn(
    match: Union[Callable, str, re.Pattern]
) -> Callable[[str], bool]:
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


def find_value(
    values,
    match: Union[Callable, str, re.Pattern],
    only_rows: Optional[List[int]] = None,
    only_cols: Optional[List[int]] = None,
) -> List[Tuple[int, int]]:
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


def cell_is_empty(cell: str):
    return _is_nan(cell)


def row_is_empty(row: List[str]) -> bool:
    """check that a row of string values is empty"""
    return all(cell_is_empty(c) for c in row)


def remove_empties_from_row(values: List[List[str]]) -> List[List[str]]:
    """Remove empty rows"""
    return [row for row in values if not row_is_empty(row)]


def is_square(values: List[List[T]]) -> bool:
    """Check that 2D list is square"""
    row_sizes = [len(r) for r in values]
    return len(set(row_sizes)) == 1


def transpose(values: List[List[T]]) -> List[List[T]]:
    """Try to 'transpose' a list of lists"""
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
    r1 = find_value(values, Keys.Basic_DNA_Parts)[0][0]
    return remove_empties_from_row(values[:r1])


def extract_basic_parts(values: List[List[Union[str, int]]]) -> List[List[CellValue]]:
    r1 = find_value(values, Keys.Basic_DNA_Parts)[0][0]
    r2 = find_value(values, Keys.Composite_DNA_Parts)[0][0]
    return remove_empties_from_row(values[r1:r2])


def extract_composite_parts(
    values: List[List[Union[str, int]]]
) -> List[List[CellValue]]:
    r1 = find_value(values, Keys.Composite_DNA_Parts)[0][0]
    return remove_empties_from_row(values[r1:])


########################
# Parse extracted values
########################


def values_to_json(
    values: List[List[Union[str, int]]],
    key_column: int,
    fill_empty_key_from_above: bool = False,
):
    data = {}
    for row in values:
        key = row[key_column].strip()
        if fill_empty_key_from_above and key == "":
            key = prev_key
        data.setdefault(key, list())
        data[key].append(row[key_column + 1 :])
        prev_key = key
    return data


# TODO: rename this method
def parse_composite_parts(values: List[List[Union[str, int]]]) -> List[Dict]:
    composite = extract_composite_parts(values)

    # locate indices of new "Collections"
    indices = find_value(composite, re.compile(Keys.Collection_Name))
    rows = [None] + [i[0] for i in indices] + [None]

    # partition the values on the indices
    partitioned = list((_partition(composite, rows)))[1:]
    parsed_json_arr = []
    for p in partitioned:
        with parsing_location(p[0][0]):
            parsed_json_arr.append(values_to_json(p, 0, True))

    # TODO: Validate the data here. Check for incorrect keys under each heading

    def _make_part_list(parsed_json):
        with parsing_location(parsed_json["Name:"][0][0]):
            names = transpose(parsed_json["Name:"])

        if 'Parts:' in parsed_json:
            # TODO: should this be covered by the jsonschema?
            # if 'Part Sequence:' in parsed_json:
            #     with parsing_location(parsed_json['Parts:'][0][0]):
            #         raise BuildRequestParsingException("Both 'Parts:' and 'Part Sequence:' cannot be "
            #                                            "defined for a composite part")
            with parsing_location(parsed_json["Parts:"][0][0]):
                parts = transpose(parsed_json["Parts:"])
            part_type = 'composite part'
        # if 'Part Sequence:' in parsed_json:
        #     with parsing_location(parsed_json['Part Sequence:'][0][0]):
        #         parts = transpose(parsed_json['Part Sequence:'])
        #     part_type = 'composite part'
        with parsing_location(parsed_json["Description:"][0][0]):
            descriptions = transpose(parsed_json["Description:"])
        collection_name = parsed_json["Collection Name:"][0][0]

        parts = [[p for p in row if not _is_nan(p)] for row in parts]

        part_list = []
        for _name, _description, _parts in list(zip(names, descriptions, parts)):
            part_list.append(
                {
                    "name": _name[0],
                    "collection": collection_name,
                    "partType": part_type,
                    "description": _description[0],
                    "parts": _parts,
                }
            )
        return part_list

    part_list = []
    for p in parsed_json_arr:
        part_list += _make_part_list(p)

    # remove empty parts (again)
    part_list_old = part_list[:]
    part_list = []
    for part in part_list_old:
        if (
            _is_nan(part['name']) and
            _is_nan(part['description']) and
                part["parts"] in [[], ["nan"], ["none"], ["None"], [""]]
        ):
            continue
        else:
            part_list.append(part)

    return part_list


def _get_and_strip(d: Dict[str, str], k: str, default: T) -> Union[T, str]:
    """Get value from dictionary and strip. Else, return default"""
    if k in d:
        if not isinstance(d[k], str):
            return default
        else:
            cv = CellValue(d[k].strip())
            cv.add_rc(d[k].row, d[k].col)
            return cv


def parse_basic_parts(values: List[List[Union[str, int]]]) -> List[Dict]:
    basic = extract_basic_parts(values)

    # remove empty rows
    basic = remove_empties_from_row(basic)

    part_json = pd.DataFrame(basic[2:], columns=basic[1]).T.to_dict()

    part_json_list = []
    for data in part_json.values():
        new_data = {
            "name": _get_and_strip(data, Keys.Part_Name, None),
            "collection": "basic DNA parts",
            "partType": "basic part",
            "description": _get_and_strip(data, Keys.Description, None),
            "sequence": _get_and_strip(data, Keys.Sequence, None),
        }

        length = data[Keys.Length]
        if length and not _is_nan(length):
            try:
                new_data["length"] = int(length)
            except ValueError:
                pass

        source = data[Keys.Source]
        if source and not _is_nan(source):
            new_data["source"] = source.strip()

        role = data[Keys.Role]
        if role and not _is_nan(role):
            new_data["roles"] = [role.strip()]

        part_json_list.append(new_data)
    return part_json_list


def parse_parts(values: List[List[Union[str, int]]]) -> List[Dict]:
    basic_parts = parse_basic_parts(values)
    composite_parts = parse_composite_parts(values)
    all_parts = basic_parts + composite_parts
    validate_part_list(all_parts)
    return all_parts
