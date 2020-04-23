import pytest

from aqbt.aquarium.parsers.pgrr_plasmid import parse_name


TEST_NAMES = [
    "pMOD-LTR1-Nat-pGRR-W5W8-URGR-F1",
    "pMOD6-pGRR-F1-yeGFP",
    "PMOD6-PGRR-F1-yeGFP",
    "pMOD8-pGRR-W5W8-iRGR-W36",
    "pMOD8-pGRR-W8-iRGR-W36",
]


@pytest.mark.parametrize("name", TEST_NAMES)
def test_parse(name):
    seq = parse_name(name)

    assert seq is not None


@pytest.mark.parametrize("name", ["pMOD8A-RGR-W36-tSynth7"])
def test_parse_fail(name):
    seq = parse_name(name)
    assert not seq


# @pytest.mark.parametrize("name", TEST_NAMES)
# def test_parse_and_upload_to_benchling(name):
#     seq = parse_name(name)
#
#     assert seq is not None
#     registry = sessions.klregistry
#     api = registry.benchling
#
#     folder = api.Folder.find_by_name("APITrash")
#     dna = convert(
#         seq, to="DNASequence", benchling_session=api, benchling_folder_id=folder.id
#     )
#
#     dna.merge(on=["name", "folder_id"])
#     print(dna.web_url)


# def test_parse_and_register():
#
#     session = sessions.aqproduction
#     with session.with_cache(timeout=60) as sess:
#         plasmid_type = sess.SampleType.find_by_name("Plasmid")
#         plasmids = sess.Sample.where({"sample_type_id": plasmid_type.id})
#         for p in plasmids:
#             seq = parse_name(p.name)
