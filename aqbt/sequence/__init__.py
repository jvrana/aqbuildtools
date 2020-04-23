"""sequence.

Methods for DNA/RNA sequence manipulation.
"""
from typing import List
from typing import Tuple

import numpy as np
from primer3plus.utils import anneal as p3panneal

BASES = "AGCT"


def random_sequence(size: int, letters=BASES) -> str:
    """Generate a random sequence string."""
    NP_BASES = np.array(list(letters))
    indices = np.random.randint(0, 4, size=size)
    bases = NP_BASES[indices]
    return "".join(bases)


def c(seq: str):
    """Return complement a dna string."""
    d = dict(zip("AGCTagct", "TCGAtcga"))
    return "".join([d[x] for x in seq])


def rc(seq: str):
    """Return a reverse complement a dna string."""
    return c(seq)[::-1]


def dna_like(seqstr: str, letters: str = "AGTCagctnNuU", min_length: int = 10):
    if seqstr is None:
        return False
    if len(seqstr) <= min_length:
        return False
    for s in seqstr:
        if s not in letters:
            return False
    return True


def anneal(
    template: str, primers: List[str], ignore_case: bool = True, n_bases: int = 10
) -> Tuple[List[dict], List[dict]]:
    """Anneal primers to a template.

    Example output.

    .. code-block:: python

        [{
            'name': None,
            'anneal': 'GTTTGCCCGGATCATCAACG',
            'overhang': '',
            'primer': 'GTTTGCCCGGATCATCAACG',
            'primer_index': 1,
            'start': 219,
            'length': 20,
            'top_strand_slice': (200, 220),
            'strand': -1
        }]

    :param template: template sequence as a string
    :param primers: list of primer sequences as a list of strings
    :param ignore_case: whether to ignore case of the primers
    :param n_bases: min num of bases to consider for annealing.
    :return: list of dictionary results
    """
    assert isinstance(template, str)
    for p in primers:
        assert isinstance(p, str)

    seq_to_index = {p.upper(): i for i, p in enumerate(primers)}
    results = p3panneal(template, primers, ignore_case=ignore_case, n_bases=n_bases)
    for r in results[0] + results[1]:
        r["primer_index"] = seq_to_index[r["primer"].upper()]
    return results
