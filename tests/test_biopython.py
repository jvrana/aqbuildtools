from aqbt import biopython
from aqbt.sequence import rc
import pytest


def test_new_sequence():
    seq = biopython.random_record(100)
    seq2 = biopython.random_record(100)
    assert len(seq) == 100
    assert len(seq) == 100
    assert sorted(list(set(seq.seq))) == ["A", "C", "G", "T"]
    assert str(seq.seq) != str(seq2.seq)


def test_annotate_sequence():
    seq = biopython.random_record(100)
    assert len(seq.features) == 1
    biopython.annotate(seq, "mylabel")
    assert len(seq.features) == 2
    for feature in seq.features:
        print(feature)

    assert int(seq.features[-1].location.start) == 0
    assert int(seq.features[-1].location.strand) == 1
    assert int(seq.features[-1].location.end) == 100
    assert seq.features[-1].qualifiers["label"][0] == "mylabel"


def test_rc():
    seq = biopython.random_record(100)
    seqrc = seq.reverse_complement()
    assert str(seqrc.seq) == rc(str(seq.seq))


def test_pcr_amplify():
    template = biopython.random_record(1000, name="template")
    primer1 = template[0:20]
    primer2 = template[200:220].reverse_complement()

    products = biopython.pcr_amplify(primer1, primer2, template, False)
    assert products
    assert len(products) == 1
    assert len(products[0][0].features) == 4


def test_no_pcr_amplify():
    template = biopython.random_record(1000, name="template")
    primer1 = template[0:20]
    primer2 = template[200:220]

    products = biopython.pcr_amplify(primer1, primer2, template, False)
    assert len(products) == 0


def test_no_pcr_amplify_over_origin():
    template = biopython.random_record(1000, name="template")
    primer1 = template[800:820]
    primer2 = template[200:220].reverse_complement()

    products = biopython.pcr_amplify(primer1, primer2, template, False)
    assert len(products) == 0


def test_pcr_amplify_over_origin():
    template = biopython.random_record(1000, name="template")
    primer1 = template[800:820]
    primer2 = template[200:220].reverse_complement()

    products = biopython.pcr_amplify(primer1, primer2, template, True)
    assert len(products) == 1


def test_annotate_cyclic():
    record = biopython.random_record(1000)
    assert len(record.features) == 1
    biopython.annotate(record, "name", 500, 200, cyclic=True)
    assert len(record.features) == 2
    assert list(record.features[-1].location) == list(range(500, 1000)) + list(
        range(200)
    )


def test_annotate_cyclic_with_non_cyclic_raises_value_err():
    record = biopython.random_record(1000)
    assert len(record.features) == 1
    with pytest.raises(ValueError):
        biopython.annotate(record, "name", 500, 200, cyclic=False)


def test_make_cyclic_assemblies():

    expected_sequence = biopython.random_record(4000, name="total")

    r1 = expected_sequence[:1000]
    r2 = expected_sequence[1000 - 30 : 2000]
    r3 = expected_sequence[2000 - 29 : 3000]
    r4 = expected_sequence[3000 - 28 : 4000] + r1[:40]

    biopython.annotate(r1, "frag1")
    biopython.annotate(r2, "frag2")
    biopython.annotate(r3, "frag3")
    biopython.annotate(r4, "frag4")

    records = [r1, r2, r3, r4]

    assemblies = biopython.make_cyclic_assemblies(records)

    # check num assemblies
    assert len(assemblies) == 1

    # check topology
    assert assemblies[0].annotations["topology"] == "circular"

    # check sequence
    assert len(assemblies[0].seq) == len(expected_sequence.seq)
    assert str(assemblies[0].seq).upper() == str(expected_sequence.seq).upper()

    # check number of features
    original_features = []
    for r in records:
        original_features += r.features

    feature_list = {}
    for f in assemblies[0].features:
        print(f)
        feature_list[f.qualifiers["label"][0]] = f

    for f in original_features:
        assert f.qualifiers["label"][0] in feature_list


def test_make_cyclic_assemblies_fail():

    expected_sequence = biopython.random_record(4000, name="total")

    r1 = expected_sequence[:1000]
    r2 = expected_sequence[970:2000]
    r3 = expected_sequence[1970:3000]
    r4 = expected_sequence[2970:3970]

    biopython.annotate(r1, "frag1")
    biopython.annotate(r2, "frag2")
    biopython.annotate(r3, "frag3")
    biopython.annotate(r4, "frag4")

    records = [r1, r2, r3, r4]

    assemblies = biopython.make_cyclic_assemblies(records)
    assert not assemblies


def test_():
    record = biopython.random_record(1000, name="random")

    record = biopython.slice_with_features(record, slice(4, 900))

    for f in record.features:
        print(f)
