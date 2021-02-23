"""Utils.

Utility methods
"""
import random
from copy import deepcopy
from itertools import tee
from typing import Any
from typing import Dict
from typing import Generator
from typing import Iterable
from typing import List
from typing import TypeVar

T = TypeVar("T")


def random_color():
    """Make a random color."""
    random_number = random.randint(0, 16777215)
    hex_number = str(hex(random_number))[2:]
    if len(hex_number) < 6:
        hex_number = "0" * (6 - len(hex_number)) + hex_number
    return "#" + hex_number


def format_float(a, places=2):
    return "%." + str(places) + "f" % round(a, places)


def remove_none_values(data: dict):
    to_remove = []
    for k, v in data.items():
        if v is None:
            to_remove.append(k)
    for k in to_remove:
        del data[k]
    return data


def random_slices(mn: int, mx: int, total_length: int):
    """Yield random slices whose lengths sum to the provided total length.

    :param mn: minimum slice size
    :param mx: maximum slice size
    :param total_length: total length of the slices
    :return:
    """
    n = total_length
    j = 0
    while j < n:
        i = j
        j += random.randint(mn, mx)
        if j >= n:
            j = n
        yield (i, j)


def sort_cycle(arr, key=None):
    """Sort a cyclic array, maintaining order."""
    if key is None:
        arr_with_i = sorted([(x, i) for i, x in enumerate(arr)])
    else:
        arr_with_i = sorted([(key(x), i) for i, x in enumerate(arr)])
    i = arr_with_i[0][1]
    return arr[i:] + arr[:i]



def chunkify(arr: Iterable[T], chunk_size: int) -> Generator[List[T], None, None]:
    new_list = []
    for x in tee(arr, 1)[0]:
        new_list.append(x)
        if len(new_list) == chunk_size:
            yield new_list
            new_list = []
    if new_list:
        yield new_list


def merge_left(
    left_dict: List[Dict[str, Any]],
    right_dict: List[Dict[str, Any]],
    left_key: str,
    right_key: str,
    new_key: str,
    inplace: bool = True,
    default: Any = None,
    merged_only: bool = False,
) -> List[Dict[str, Any]]:
    """Merge the list of dictionaries with another list of dictionaries on the
    provided keys by matching `left_key` and `right_key` keys and creating a
    new attribute `new_key`

    .. code-block:: python
        new_key = "new_key"

        left = {
            "key1": 1
            "left_key": 3
        }

        right = {
            "key2": 2,
            "right_key": 3
        }

        merged = merge_left(left, right, "left_key", "right_key", "new_key")

        assert merged == {
            "key1": 1,
            "left_key": 3,
            "new_key": [{
                "key2": 2,
                "right_key": 3
            }]
        }

    :param left_dict: left dict to merge
    :param right_dict: right dict to merge
    :param left_key: key of the left dict to merge
    :param right_key: key of the right dict to merge
    :param new_key: new attribute
    :param inplace: if False, return copy of the dictionary
    :param default: default value to add if merge
    :param merged_only: return only the merged di
    :return: merged list of dictionaries
    """
    if not inplace:
        left_dict = deepcopy(left_dict)
    d1_key_dict = {entry[left_key]: entry for entry in left_dict}
    d2_key_dict = {entry[right_key]: entry for entry in right_dict}
    for k, v in d1_key_dict.items():
        d1_key_dict[k][new_key] = default
    if merged_only:
        merged = []
    for k, v in d2_key_dict.items():
        if k in d1_key_dict:
            x = d1_key_dict[k]
            if merged_only:
                merged.append(x)
            x.setdefault(new_key, list())
            if x[new_key] is None:
                x[new_key] = []
            x[new_key].append(v)
    if merged_only:
        return merged
    return left_dict
