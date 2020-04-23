# pgrr_plasmid.py
#
# specialized parser function for Gander 2017 style plasmid names
import functools
import json
import os
import re
from copy import deepcopy
from typing import Any
from typing import List
from typing import Union

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from pydent import AqSession

from aqbt import biopython

# import jdna
# from colorhash import ColorHash

here = os.path.abspath(os.path.dirname(__file__))
data_dir = os.path.join(here, "data")


def parts_df():
    PARTS_FILEPATH = os.path.join(data_dir, "NORGateParts.csv")

    parts_df = pd.read_csv(PARTS_FILEPATH)
    parts_df.fillna(value="None", inplace=True)
    parts_df.tail()
    return parts_df


def parts_by_name_and_cat():
    parts_dict = {}
    by_category = {}
    df = parts_df()
    for _, row in df.iterrows():
        seq = row["Sequence"]
        if seq != "None":
            part_name = row["Part"]
            seq = biopython.new_sequence(
                row["Sequence"],
                name=part_name,
                auto_annotate=True,
                cyclic=False,
                annotation_type=row["Annotation"],
            )
            parts_dict[part_name.lower()] = seq
            by_category.setdefault(row["Annotation"], {})
            by_category[row["Annotation"]][part_name.lower()] = seq
    return parts_dict, by_category


_parts_dict, _parts_by_category_dict = parts_by_name_and_cat()


class Null:
    pass


def get_part(part: Union[str, SeqRecord], default: Any = Null):
    if isinstance(part, str):
        part_instance = _parts_dict.get(part.lower(), default)
        if part_instance is Null:
            return _parts_dict[part.lower()]
        return deepcopy(part_instance)
    return part


def wsetdict():
    """Return a dictionary converting w-set to r-set notation.

    :return:
    """
    path = os.path.join(data_dir, "wset_conversion.json")
    with open(path, "r") as f:
        wset = json.load(f)
    return wset


def concat_parts(parts: List[Union[str, SeqRecord]]) -> SeqRecord:
    """Concatenates a list of Coral.DNA objects or looks for parts by name in
    the parts dictionary."""
    record = biopython.new_sequence("", name="")
    for p in parts:
        if p is None:
            continue
        if isinstance(p, str):
            p = p.strip()
            if p == "":
                continue
            record += get_part(p)
        else:
            record += p
    return record


def generate_gRNA_cassette(
    prom: Union[str, SeqRecord] = None,
    promlink: Union[str, SeqRecord] = None,
    target: Union[str, SeqRecord] = None,
    handle: Union[str, SeqRecord] = None,
    ribo5: Union[str, SeqRecord] = None,
    ribo3: Union[str, SeqRecord] = None,
    ribo5INS=(None, None),
    termlink: Union[str, SeqRecord] = None,
    term: Union[str, SeqRecord] = None,
):
    """Generates gRNA cassette DNA sequence using Coral."""
    parts = [
        prom,
        promlink,
        ribo5INS[0],
        ribo5,
        ribo5INS[1],
        target,
        handle,
        ribo3,
        termlink,
        term,
    ]
    seq = concat_parts(parts)
    return seq


# helper function for uninsulated RGR cassette
helper_rgr = functools.partial(
    generate_gRNA_cassette,
    handle="spCas9 gRNA Handle",
    ribo3="HDV Ribozyme",
    ribo5="HH Ribozyme",
    termlink="TP",
    term="tCYC1",
)


# uninsulated RGR cassette
def rgr(target):
    ins = deepcopy(target)[:6].reverse_complement()
    return helper_rgr(target=deepcopy(target), ribo5INS=(ins, None))


# insulated pADH1 gRNA cassette
padh1_irgr = functools.partial(
    generate_gRNA_cassette,
    handle="spCas9 gRNA Handle",
    ribo3="HDV Ribozyme",
    ribo5="ASBV1 Ribozyme",
    ribo5INS=("ASBV1 pADH1_INS1", "ASBV1 pADH1_INS2"),
    termlink="TP",
    term="tCYC1",
)

# insulated pGAZL4 gRNA cassette
pgalz4_irgr = functools.partial(
    generate_gRNA_cassette,
    handle="spCas9 gRNA Handle",
    ribo3="HDV Ribozyme",
    ribo5="ASBV1 Ribozyme",
    ribo5INS=("ASBV1 pGAL1_INS1", "ASBV1 pGAL1_INS2"),
    termlink="TP",
    term="tCYC1",
)

# insulated pGRR gRNA cassette
pGRR_irgr = functools.partial(
    generate_gRNA_cassette,
    handle="spCas9 gRNA Handle",
    ribo3="HDV Ribozyme",
    ribo5="ASBV1 Ribozyme",
    ribo5INS=("ASBV1 pGRR_INS1", "ASBV1 pGRR_INS2"),
    termlink="TP",
    term="tCYC1",
)


def pGRR(i: str, j: str) -> SeqRecord:
    """Generates a new pGRR promoter from inputs 'i' and 'j'."""
    i = get_part(i, default=None)
    j = get_part(j, default=None)
    if i is not None:
        i = deepcopy(i)
    if j is not None:
        j = deepcopy(j)
        j = j.reverse_complement()

    parts = ["pGRR", i, "pGRR TATA", j, "pGRR RBS"]
    return concat_parts(parts)


def pMOD(homology: str, cassette: SeqRecord, marker=None):
    """Generates new pMOD vector with a new cassette.

    Cassette should include the promoter, gene, and terminator (no
    assembly linkers).
    """

    homology_dict = {
        "URA": [["URA3 Promoter", "URA3", "tADH1"], ["URA3 3'UTR"]],
        "HIS": [["HIS3 Promoter", "HIS3", "tADH1"], ["HIS3 3'UTR"]],
        "TRP": [["TRP1 Promoter", "TRP1", "tADH1"], ["TRP1 3'UTR"]],
        "LTR1": [["LTR1 homology 1", "NATMX"], ["LTR1 homology 2"]],
        "LTR2": [["LTR2 homology 1", "BLEOMX"], ["LTR2 homology 2"]],
        "LTR3": [["LTR3 homology 1", "HYBMX"], ["LTR3 homology 2"]],
        "HO": [["HO homology 1", "KANMX"], ["HO homology 2"]],
    }

    hom1, hom2 = homology_dict[homology]
    if marker is not None:
        hom1[1] = marker

    _parts = ["PmeI", hom1, "PP2", cassette, "TS", hom2, "PmeI", "AMPR and ORI"]
    _parts = _parts[1:-2]
    parts = []
    for p in _parts:
        if isinstance(p, list):
            parts += p
        else:
            parts.append(p)
    return concat_parts(parts)


def parse_NOR_gate(name: str) -> re.Match:
    # {'marker': '8', 'prom': 'A', 'i': None, 'j': None, 'grna': 'URGR', 'k': 'W10'}
    imatch = r"(?P<i>[WF]\d+)"
    jmatch = r"(?P<j>[WF]\d+)"

    pgrr_match = r"(pGRR(-{i})?(-?{j})?)".format(i=imatch, j=jmatch)
    promoter_match = r"p\w+|{pgrr_match}|[GA]".format(pgrr_match=pgrr_match)

    marker_match = (
        r"(?P<marker_int>\d)|(?P<homology>ho|HO|\w+)-?(?P<marker>kan|KAN|Kan|\w+)"
    )

    cassette_match = r"(?P<grna>[iU]?RGR)-(?P<k>[WF]\d+)|[y]?eGFP"
    pattern = r"pMOD-?({marker})-?(?P<prom>{promoter})-(?P<cassette>{cassette})$".format(
        cassette=cassette_match, marker=marker_match, promoter=promoter_match
    )
    m = re.match(pattern, name, re.IGNORECASE)
    return m


amp_ori = SeqIO.read(os.path.join(here, "data", "amporibackbone.gb"), format="genbank")


def parse_nor_gate_name_to_sequence(name: str) -> Union[SeqRecord, None]:
    match = parse_NOR_gate(name)
    if match is None:
        return None

    parsed = match.groupdict()
    homology_dict = {6: "URA", 8: "HIS", 4: "TRP"}
    promoter_dict = {"A": "pADH1", "G": "pGPD"}
    irgr_dict = {"pADH1": padh1_irgr, "pGALZ4": pgalz4_irgr, "pGRR": pGRR_irgr}

    def get_target(x):
        return wsetdict().get(x, x)

    i = get_target(parsed["i"])
    j = get_target(parsed["j"])
    k = get_target(parsed["k"])

    if i is None:
        i = "c3"
    if j is None:
        j = "c6"

    marker_int = parsed["marker_int"]
    if marker_int is not None:
        marker_int = int(marker_int)

    homology = homology_dict.get(marker_int, parsed["homology"])
    if homology:
        try:
            homology = homology_dict[int(homology)]
        except ValueError:
            pass
    else:
        return None
    marker = parsed["marker"]
    prom = promoter_dict.get(parsed["prom"], parsed["prom"])

    grna = None
    if parsed["grna"] is not None:
        if parsed["grna"] == "RGR":
            grna = rgr(get_part(k))
        elif parsed["grna"].lower() in ["irgr", "urgr"]:
            if "pGRR".upper() in prom.upper():
                grna = irgr_dict["pGRR"](k)
            else:
                grna = irgr_dict[prom](k)

    if "pGRR".lower() in prom.lower():
        prom = pGRR(i, j)

    cassette = parsed["cassette"]
    if grna is not None:
        cassette = concat_parts([prom, "PS", grna])
    else:
        if "egfp" in cassette.lower():
            cassette = "yeGFP"
        cassette = [prom, "PS", cassette, "TP", "tCYC1"]

    if marker is not None:
        marker = marker.upper()
        if "hyg" in marker.lower():
            marker = "HYGMX"
        if "nat" in marker.lower():
            marker = "NATMX"
        if "kan" in marker.lower():
            marker = "KANMX"
        if "bleo" in marker.lower():
            marker = "BLEOMX"

    record = pMOD(homology.upper(), cassette, marker=marker)
    biopython.make_cyclic([record])
    record.id = name
    record.annotations["name"] = name
    record.name = name

    record = record + amp_ori
    biopython.remove_duplicate_features(record)

    return record


# alias for main function


def parse_name(name: str) -> Union[None, SeqRecord]:
    try:
        return parse_nor_gate_name_to_sequence(name)
    except KeyError:
        return None
