from .cell_value import CellValue
from typing import Union


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
            exception = BuildRequestParsingException(
                "{} near ({}, {}).\n{}".format(exc_type, self.row, self.col, exc_value)
            )
            exception.row = self.row
            exception.col = self.col
            raise exception


class EmptyContext(object):
    def __enter__(self):
        pass

    def __exit__(self, exc_type, exc_val, exc_tb):
        pass


def parsing_location(value: Union[str, CellValue]):
    if isinstance(value, CellValue):
        return LocationContext(*value.rc())
    else:
        return EmptyContext()
