import pytest
from aqbt import biopython
from os.path import join, abspath, dirname
from aqbt.aquarium.genome_builder import integration, mating


@pytest.fixture(scope='function')
def cenpk(fixtures_path):
    return biopython.from_gffs([join(fixtures_path, 'cenpk.gff')])


@pytest.mark.parametrize(
    "sample_id",
    [
        27351,
        27673,
        27333,
        27683,
        33174,
        33175,
        33176,
        33177,
        33178,
        33179,
        33180,
        33181,
    ],
)
def test_integrations(sample_id, cenpk, registry, tmp_path):
    records = cenpk

    d = {22800: {"genome": records[:]}, 22801: {"genome": records[:]}}
    session = registry.session

    genome_dict = integration(
        registry, session.Sample.find(sample_id), d, "PmeI"
    )

    for sid, genome in genome_dict.items():
        name = "UWBF_{}.gff".format(sid)
        biopython.to_gff(genome["genome"], join(tmp_path, name))


def test_mating(cenpk):
    mata = cenpk[:]
    alpha = cenpk[1:]
    records = mating(mata, alpha)
    assert len(records) == 2*len(mata) - 1

# @pytest.mark.parametrize(
#     "sample_id",
#     [
#         27351,
#         27673,
#         27333,
#         27683,
#         33174,
#         33175,
#         33176,
#         33177,
#         33178,
#         33179,
#         33180,
#         33181,
#     ],
# )
# def test_aq_to_gff(sample_id):
#     records = list(GFF.parse([join(here, "cenpk.gff")]))
#
#     d = {22800: {"genome": records[:]}, 22801: {"genome": records[:]}}
#
#     registry = sessions.klregistry
#     session = registry.session
#
#     path, genome_dict, trace = aq_to_gff(
#         registry, session.Sample.find(sample_id), d, "PmeI", do_save=False
#     )
#
#
# def test_aq_to_gff_example():
#     records = list(GFF.parse([join(here, "cenpk.gff")]))
#
#     d = {22800: {"genome": records[:]}, 22801: {"genome": records[:]}}
#
#     registry = sessions.klregistry
#     session = registry.session
#
#     sample = session.Sample.find(34018)
#
#     path, genome_dict, trace = aq_to_gff(
#         registry, sample, d, "PmeI", do_save=True, out_dir=abspath(dirname(__file__))
#     )
#
#     for t in trace:
#         for i in t["integrations"]:
#             s = i["left_flank"][-1000:] + i["cassette"] + i["right_flank"][:1000]
#             from Bio.Seq import Alphabet
#
#             s.seq.alphabet = Alphabet.generic_dna
#             long_name = "{}__{}".format(sample.id, sample.name)
#             short_name = "UWBF_{}".format(sample.id)
#             s.name = short_name
#             s.id = short_name
#             SeqIO.write(
#                 [s],
#                 handle=join(abspath(dirname(__file__)), "{}.gb".format(long_name)),
#                 format="genbank",
#             )
#             print(
#                 {
#                     "left_flank": i["left_flank"][-1000:],
#                     "cassette": i["cassette"],
#                     "right_flank": i["right_flank"][:1000],
#                 }
#             )
#
#
# @pytest.mark.parametrize("sample_id", [27674])
# def test_aq_to_gff_diploid(sample_id):
#     records = list(GFF.parse([join(here, "cenpk.gff")]))
#
#     d = {22800: {"genome": records[:]}, 22801: {"genome": records[:]}}
#
#     registry = sessions.klregistry
#     session = registry.session
#
#     path, genome_dict = aq_to_gff(
#         registry, session.Sample.find(sample_id), d, "PmeI", do_save=False
#     )
