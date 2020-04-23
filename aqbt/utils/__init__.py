"""Utils.

Utilities for lobio.
"""
import random
from itertools import tee
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
    counter = 0
    for x in tee(arr, 1)[0]:
        if counter >= chunk_size:
            yield new_list
            new_list = []
            counter = 0
        else:
            new_list.append(x)
            counter += 1
    yield new_list
