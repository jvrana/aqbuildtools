import hashlib
import os
from glob import glob
from os.path import abspath
from os.path import dirname
from os.path import join

import pytest
from dasi import cost

##############################
# Cost Fixtures
##############################
here = dirname(abspath(__file__))


@pytest.fixture(scope="session")
def cost_filepath():
    return join(here, "span_cost.b")


@pytest.fixture(scope="session")
def cost_checksum_filepath():
    return join(here, "span_cost_checksum.txt")


def hashfiles(files, hash="sha1", hashfunc=None):
    if not hashfunc:

        def hashfunc(string):
            return getattr(hashlib, hash)(string.encode("utf-8")).hexdigest()

    contents = ""
    sorted_path = sorted(files)
    for file in sorted_path:
        with open(file, "r") as f:
            contents += hashfunc(f.read())
    return hashfunc(contents)


def cost_module_checksum():
    cost_dir = os.path.dirname(cost.__file__)
    cost_files = sorted(
        glob(os.path.join(cost_dir, "*.py")) + glob(os.path.join(cost_dir, "*.json"))
    )
    return hashfiles(cost_files)


def cached(path, save_func, load_func, checksum_path, logger=None):
    different_checksum = True
    checksum = cost_module_checksum()
    if os.path.isfile(checksum_path):
        with open(checksum_path, "r") as f:
            stored_checksum = f.read().strip()
            if stored_checksum == checksum:
                different_checksum = False
            if logger:
                logger.debug("Stored checksum: {}".format(stored_checksum))
    if logger:
        logger.debug("Checksum: {}".format(checksum))
    if different_checksum or not os.path.isfile(path):
        if logger:
            logger.debug("Using default params")
        model = save_func(path)
        with open(checksum_path, "w") as f:
            f.write(checksum)
        stat = os.stat(path)
        if logger:
            logger.debug("Wrote {} bytes".format(stat.st_size))
    else:
        if logger:
            logger.debug("Loading {}".format(path))
        stat = os.stat(path)
        if logger:
            logger.debug("Loaded {} bytes".format(stat.st_size))
        model = load_func(path)
    return model


@pytest.fixture(scope="session")
def cached_span_cost(cost_filepath, cost_checksum_filepath):
    """This will check the checksum of the cost module against the last
    checksum. If checksums are the same, the span cost will be loaded. Else,
    span_cost will be created from default parameters and saved with the cost
    module's checksum.

    :param cost_filepath: path of the span_cost
    :param cost_checksum_filepath: path of the checksum
    :return: SpanCost
    """

    def load_span_cost(path):
        span_cost = cost.SpanCost.load(path)
        return span_cost

    def save_span_cost(path):
        span_cost = cost.SpanCost.open()
        span_cost.dump(path)
        return span_cost

    return cached(
        cost_filepath,
        load_func=load_span_cost,
        save_func=save_span_cost,
        checksum_path=cost_checksum_filepath,
    )
