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