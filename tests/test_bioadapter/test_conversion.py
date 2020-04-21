from copy import deepcopy

import dictdiffer
import networkx as nx
import pytest
from benchlingapi import Session
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pydent import AqSession
from tqdm import tqdm

from aqbt.bioadapter import bioadapter


benchling = Session("sk_op57DfvaGDESTkPysNb57THFy5GO1")
session = AqSession("vrana", "Mountain5", "http://52.27.43.242/")


def edge_iter(src, as_str=False):
    shortest_paths = nx.shortest_path(bioadapter.g, source=src)
    assert shortest_paths
    for k, v in shortest_paths.items():
        if v[0] == v[-1]:
            continue
        if as_str:
            yield "{} {}".format(v[0], v[-1])
        else:
            yield v[0], v[-1]


@pytest.fixture(scope="module")
def benchling_dna_list():
    return benchling.DNASequence.last(50)


@pytest.fixture(scope="module")
def benchling_dna(benchling_dna_list):
    return benchling_dna_list[0]


@pytest.fixture(scope="module")
def aquarium_sample_list():
    st = session.SampleType.find_by_name("Plasmid")
    return session.Sample.where({"sample_type_id": st.id, "user_id": 66})


@pytest.fixture(scope="module")
def aquarium_sample(aquarium_sample_list):
    return aquarium_sample_list[0]


def parametrize_paths(src):
    edges = list(edge_iter(src))
    # ids = [str(x[0]) + " -> " + str(x[1]) for x in edges]
    ids = [str(bioadapter.path(x[0], x[1])) for x in edges]
    return pytest.mark.parametrize("edge", edges, ids=ids)


@parametrize_paths("benchling_dna_json")
@pytest.mark.parametrize("unconvert", [False, True], ids=["", "reversed"])
def test_benchling_dna_json_conversion(edge, unconvert, benchling_dna):
    data = benchling_dna.dump()
    converted = bioadapter.convert(
        data,
        edge[0],
        edge[-1],
        benchling_session=benchling,
        path="dna.gb",
        format="fasta",
    )
    if unconvert:
        unconverted = bioadapter.convert(
            converted, edge[-1], edge[0], benchling_session=benchling, format="fasta"
        )
        difference = dictdiffer.diff(data, unconverted)
        for line in difference:
            print(line)
        if edge[-1] in ["genbank", "file", "fasta"]:
            assert type(unconverted) is type(data)
        else:
            assert unconverted == data


@parametrize_paths("list(benchling_dna_json)")
@pytest.mark.parametrize("unconvert", [False, True], ids=["", "reversed"])
def test_benchling_dna_list_json_conversion(edge, unconvert, benchling_dna_list):
    data = [d.dump() for d in benchling_dna_list]
    converted = bioadapter.convert(
        data,
        edge[0],
        edge[-1],
        benchling_session=benchling,
        path="dna.gb",
        format="fasta",
    )
    print(converted)
    if unconvert:
        if edge[-1] in ["list(genbank)", "list(file)", "list(fasta)"]:
            pytest.skip("Saving files loses important information.")
        unconverted = bioadapter.convert(
            converted, edge[-1], edge[0], benchling_session=benchling, format="genbank"
        )
        difference = dictdiffer.diff(data, unconverted)
        for line in difference:
            print(line)
        assert unconverted == data


@pytest.mark.parametrize("pbar", [True, False, tqdm])
def test_pbar(pbar, benchling_dna_list):
    data = [d.dump() for d in benchling_dna_list]
    converted = bioadapter.convert(
        data,
        "list(benchling_dna_json)",
        "list(SeqRecord)",
        benchling_session=benchling,
        path="dna.gb",
        format="fasta",
        pbar=pbar,
    )
    assert converted


@parametrize_paths("DNASequence")
def test_no_raise(benchling_dna, edge):
    """no raise should return None."""
    dna = deepcopy(benchling_dna)
    del dna.__dict__["bases"]
    if edge[-1] == "benchling_dna_json":
        pytest.skip()
    converted = bioadapter.convert(dna, edge[0], edge[-1], no_raise=True)
    assert converted is None


# @parametrize_paths("list(DNASequence)")
def test_no_raise_list(benchling_dna_list):
    """no raise should return None."""
    dna = deepcopy(benchling_dna_list)
    del dna[0].__dict__["bases"]
    converted = bioadapter.convert(dna, to="SeqRecord", no_raise=True)
    assert len(converted) == len(benchling_dna_list)
    assert converted[0] is None
    for c in converted[1:]:
        assert c is not None


@parametrize_paths("benchling_annotation_json")
@pytest.mark.parametrize("unconvert", [False, True], ids=["", "reversed"])
def test_benchling_feature_json_conversion(edge, unconvert, benchling_dna):
    data = benchling_dna.annotations[0].dump()
    converted = bioadapter.convert(
        data,
        edge[0],
        edge[-1],
        length=len(benchling_dna.bases),
        benchling_session=benchling,
    )
    if unconvert:
        unconverted = bioadapter.convert(
            converted,
            edge[-1],
            edge[0],
            length=len(benchling_dna.bases),
            benchling_session=benchling,
        )
        difference = dictdiffer.diff(data, unconverted, ignore=["type"])
        for line in difference:
            print(line)
        assert unconverted == data


@parametrize_paths("Sample")
def test_aquarium_sample_conversion(edge, aquarium_sample):
    converted = bioadapter.convert(
        aquarium_sample,
        edge[0],
        edge[-1],
        benchling_session=benchling,
        path="dna.gb",
        format="fasta",
    )


@parametrize_paths("list(Sample)")
def test_aquarium_sample_list_conversion(edge, aquarium_sample_list):
    converted = bioadapter.convert(
        aquarium_sample_list[:10],
        edge[0],
        edge[-1],
        benchling_session=benchling,
        path="dna.gb",
        format="fasta",
        no_raise=True,
        pbar=True,
    )
    assert converted


@parametrize_paths("SeqRecord")
def test_seqrecord_conversion(edge):
    record = SeqRecord(Seq("AGGTAGGATTA", generic_dna), annotations={"folderId": None})
    bioadapter.convert(
        record,
        edge[0],
        edge[-1],
        benchling_session=benchling,
        path="dna.gb",
        format="genbank",
    )


def test_fasta_to_genbank():
    record = SeqRecord(Seq("AGGTAGGATTA", generic_dna), annotations={"folderId": None})
    bioadapter.convert(record, to="fasta", path="dna.fasta")
    bioadapter.convert("dna.fasta", frm="fasta", to="genbank", path="dna.gb")


def test_genbank_to_fasta():
    record = SeqRecord(Seq("AGGTAGGATTA", generic_dna), annotations={"folderId": None})
    bioadapter.convert(record, to="genbank", path="dna.gb")
    bioadapter.convert("dna.gb", frm="genbank", to="fasta", path="dna.fasta")


@parametrize_paths("genbank")
def test_fasta_conversion(edge):
    record = SeqRecord(Seq("AGGTAGGATTA", generic_dna), annotations={"folderId": None})
    path = bioadapter.convert(record, to="genbank", path="dna.gb")
    bioadapter.convert(
        path,
        edge[0],
        edge[-1],
        benchling_session=benchling,
        path="dna.gb",
        format="genbank",
        benchling_folder_id="asdfl34",
    )


@parametrize_paths("list(SeqRecord)")
def test_seqrecord_list_conversion(edge):
    record1 = SeqRecord(Seq("AGGTAGGATTA", generic_dna), annotations={"folderId": None})
    record2 = SeqRecord(
        Seq("AGasdfasdGTAGGATTA", generic_dna), annotations={"folderId": None}
    )
    record3 = SeqRecord(
        Seq("AGGTdfadsfAGGATTA", generic_dna), annotations={"folderId": None}
    )
    records = [record1, record2, record3]
    bioadapter.convert(
        records,
        edge[0],
        edge[-1],
        benchling_session=benchling,
        path="dna.gb",
        format="genbank",
    )
