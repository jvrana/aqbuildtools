"""
bioadapter

Provides methods and classes for biological data type conversion

e.g. genbank -> benchling -> fasta -> ... etc.

Registers several functions and uses the expected data type io to
create a graph to find the optimal conversion that minimizes data
loss.
"""

from importlib import import_module as _import_module
from .registry import convert
from .registry import bioadapter


_import_module("aqbt.bioadapter.conversion", "bioadapter")