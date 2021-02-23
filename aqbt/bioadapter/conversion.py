from copy import deepcopy
from typing import List

from BCBio import GFF
from benchlingapi.models import Annotation
from benchlingapi.models import DNASequence
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation
from Bio.SeqFeature import FeatureLocation
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord
from pydent.models import Sample

from ._bioadapter import bioadapter

################################################
# BioPython Adapter
################################################

benchling_annotation_json = "benchling_annotation_json"
benchling_dna_json = "benchling_dna_json"


@bioadapter(benchling_annotation_json, SeqFeature, many=True)
def json_to_feature(data: dict, length: int) -> SeqFeature:
    start, end, strand = data["start"], data["end"], data["strand"]
    if end == 0:
        end = length

    if end < start:
        location = CompoundLocation(
            [
                FeatureLocation(start, length, strand=strand),
                FeatureLocation(0, end, strand=strand),
            ]
        )
    else:
        location = FeatureLocation(start, end, strand=strand)
    feature = SeqFeature(
        location=location,
        type=data.get("type", "misc_feature"),
        qualifiers={
            "label": [data["name"]],
            "ApEinfo_fwdcolor": [data["color"]],
            "ApEinfo_revcolor": [data["color"]],
        },
    )
    if feature.type.strip() == "":
        feature.type = "misc_feature"
    return feature


@bioadapter(SeqFeature, benchling_annotation_json, many=True)
def feature_to_json(feature: SeqFeature, length: int) -> dict:
    location = feature.location
    if isinstance(location, CompoundLocation):
        start = location.parts[0].start
        end = location.parts[-1].end
        strand = location.parts[0].strand
        assert location.parts[0].end == length
    else:
        start = int(feature.location.start)
        end = int(feature.location.end)
        strand = feature.location.strand
    qualifier_label = feature.qualifiers["label"]
    if isinstance(qualifier_label, list):
        name = qualifier_label[0]
    elif isinstance(qualifier_label, str):
        name = qualifier_label
    elif qualifier_label is None:
        name = "misc_feature"
    else:
        raise ValueError(
            "Type '{}' {} is not a valid feature name".format(
                type(qualifier_label), qualifier_label
            )
        )
    color = feature.qualifiers["ApEinfo_fwdcolor"][0]
    if end >= length:
        end = 0
    if start >= length:
        start = 0
    return {
        "start": start,
        "end": end,
        "strand": strand,
        "color": color,
        "name": name,
        "type": feature.type,
    }


def _get_benchling_dna_json_annotations(data: dict) -> dict:
    annotations = {
        "fields": None,
        "creator": None,
        "folderId": None,
        "webUrl": None,
        "customFields": None,
        "translations": None,
        "schema": None,
        "aliases": None,
        "registryId": None,
        "entityRegistryId": None,
        "schemaId": None,
    }
    for k in annotations:
        if k in data:
            annotations[k] = deepcopy(data[k])
    return annotations


@bioadapter("benchling_dna_json", SeqRecord, many=True)
def json_to_seqrecord(data: dict) -> SeqRecord:
    features = [json_to_feature(a, data["length"]) for a in data["annotations"]]
    topology = "linear"
    if data["isCircular"]:
        topology = "circular"
    annotations = _get_benchling_dna_json_annotations(data)
    annotations["topology"] = topology
    return SeqRecord(
        Seq(data["bases"]),
        name=data["name"],
        id=data.get("id", None),
        annotations=annotations,
        features=features,
    )


@bioadapter(SeqRecord, benchling_dna_json, many=True)
def seqrecord_to_json(record: SeqRecord) -> dict:
    feature_annotations = [feature_to_json(f, len(record.seq)) for f in record.features]
    data = {
        "bases": str(record.seq),
        "length": len(record.seq),
        "name": record.name,
        "id": record.id,
        "annotations": feature_annotations,
    }
    annotations = deepcopy(record.annotations)

    if annotations.get("topology", "linear") == "circular":
        annotations["isCircular"] = True
    else:
        annotations["isCircular"] = False
    if "topology" in annotations:
        del annotations["topology"]
    data.update(annotations)
    return data


@bioadapter(SeqRecord, "genbank", weight=10)
def seqrecord_to_genbank(seqrecord: SeqRecord, path: str) -> str:
    seqrecord = deepcopy(seqrecord)
    seqrecord.name = seqrecord.name.replace(" ", "")
    SeqIO.write(seqrecord, path, format="genbank")
    return path


@bioadapter("genbank", SeqRecord, weight=10)
def genbank_to_seqrecord(path: str) -> SeqRecord:
    return SeqIO.read(path, format="genbank")


@bioadapter("list(SeqRecord)", "list(genbank)", weight=10)
def seqrecords_to_genbank(seqrecordlist: List[SeqRecord], path: str) -> str:
    records = []
    for seq in seqrecordlist:
        seq = deepcopy(seq)
        seq.name = seq.name.replace(" ", "")
        records.append(seq)
    SeqIO.write(records, path, format="genbank")
    return path


@bioadapter("list(genbank)", "list(SeqRecord)", weight=10)
def genbank_to_seqrecords(path: str) -> List[SeqRecord]:
    seqs = []
    for seq in SeqIO.parse(path, format="genbank"):
        seqs.append(seq)
    return seqs


@bioadapter("SeqRecord", "fasta", weight=15)
def seqrecord_to_fasta(record: SeqRecord, path: str) -> str:
    SeqIO.write([record], path, format="fasta")
    return path


@bioadapter("fasta", "SeqRecord", weight=15)
def fasta_to_seqrecord(path: str) -> SeqRecord:
    record = SeqIO.read(path, format="fasta")
    return record


@bioadapter("list(fasta)", "list(SeqRecord)", weight=15)
def fasta_to_seqrecords(path: str) -> List[SeqRecord]:
    seqs = []
    for seq in SeqIO.parse(path, format="fasta"):
        seqs.append(seq)
    return seqs


@bioadapter("list(SeqRecord)", "list(fasta)", weight=15)
def seqrecords_to_fasta(records: List[SeqRecord], path: str) -> str:
    SeqIO.write(records, path, format="fasta")
    return path


@bioadapter("SeqRecord", "file", weight=30)
def seqrecord_to_file(record: SeqRecord, path: str, format: str) -> str:
    SeqIO.write([record], path, format=format)
    return path


@bioadapter("file", "SeqRecord", weight=30)
def file_to_seqrecord(path: str, format: str) -> SeqRecord:
    return SeqIO.read(path, format=format)


@bioadapter("list(file)", "list(SeqRecord)", weight=30)
def files_to_seqrecords(path: str, format: str) -> SeqRecord:
    seqs = []
    for seq in SeqIO.parse(path, format=format):
        seqs.append(seq)
    return seqs


@bioadapter("list(SeqRecord)", "list(file)", weight=10)
def seqrecords_to_file(seqrecordlist: List[SeqRecord], path: str, format: str) -> str:
    SeqIO.write(seqrecordlist, path, format=format)
    return path


# @bioadapter("list(SeqRecord)", "gff", weight=10, many=False)
# def seqrecords_to_gff(seqrecordlist: List[SeqRecord], path: str):
#     with open(path, 'w') as f:
#         GFF.write(seqrecordlist, f, include_fasta=True)
#     return str(path)
#
# @bioadapter("list(gff)", "list(SeqRecord)", weight=10, many=False)
# def gff_to_seqrecords(gffs: List[str]) -> List[SeqRecord]:
#     return list(GFF.parse(gffs))


@bioadapter("gff", "list(SeqRecord)", weight=10, many=False)
def gff_to_seqrecords(gff: str) -> List[SeqRecord]:
    return list(GFF.parse([gff]))


################################################
# Benchling Adapter
################################################


@bioadapter(benchling_annotation_json, Annotation, many=True)
def json_to_annotation(data: dict) -> Annotation:
    return Annotation(**data)


@bioadapter(Annotation, benchling_annotation_json, many=True)
def annotation_to_json(annotation: Annotation) -> dict:
    return annotation.dump()


@bioadapter(DNASequence, benchling_dna_json, many=True)
def benchling_to_json(dna: DNASequence) -> dict:
    return dna.dump()


@bioadapter(benchling_dna_json, DNASequence, many=True)
def json_to_benchling(data: dict, benchling_session, benchling_folder_id=None):
    if benchling_folder_id:
        data = deepcopy(data)
        data["folder_id"] = benchling_folder_id
    return benchling_session.DNASequence.load(data)


@bioadapter("benchling_sharelink", DNASequence, many=True)
def sharelink_to_benchling(sharelink: str, benchling_session) -> DNASequence:
    if not sharelink.strip():
        return None
    dna = benchling_session.DNASequence.from_share_link(sharelink.strip())
    return dna


################################################
# jDNA Adapter
################################################

# @bioadapter(benchling_dna_json, jDNASequence, weight=10)
# def json_to_jdna(data: dict) -> jDNASequence:
#
#     metadata = _get_benchling_dna_json_annotations(data)
#     seq = jDNASequence(data['bases'],
#                  name=data['name'],
#                  description=data['description'],
#                  cyclic=data["isCircular"],
#                 metadata=metadata)
#     for a in data['annotations']:
#         feature = jDNAFeature(name=a['name'], type=a['type'], strand=a['strand'], color=a['color'])
#         seq.add_feature(a['start'], a['end'], feature)
#     return seq
#
# @bioadapter(jDNASequence, benchling_dna_json, weight=10)
# def jdna_to_json(data: jDNASequence) -> dict:
#     pass
# {
#         "start": start,
#         "end": end+1,
#         "strand": strand,
#         "color": color,
#         "name": name,
#         "type": feature.type,
#     }


################################################
# Aquarium Adapter
################################################


@bioadapter(Sample, "benchling_sharelink")
def aquarium_to_sharelink(sample: Sample) -> DNASequence:
    for fv in sample.field_values:
        if fv.name == "Sequence":
            weburl = fv.value
            return weburl


@bioadapter("list(Sample)", "list(benchling_sharelink)")
def aquariums_to_sharelinks(samples: List[Sample]) -> List[str]:
    with samples[0].session.with_cache(using_models=True, timeout=60) as sess:
        sess.browser.get(samples, "field_values")
        share_links = []
        for sample in samples:
            prop = None
            for fv in sample.field_values:
                if fv.name == "Sequence":
                    prop = fv
                break
            if prop:
                share_links.append(fv.value.strip())
            else:
                share_links.append(None)
    return share_links


# @registry(benchling_dna_json, Sample)
# def json_to_aquarium(json: dict, aq_session: AqSession) -> Sample:
#     g = re.search(r'aqUWBF(\d+)', json['entity_registry_id'])
#     aqid = int(g.group(1))
#     return aq_session.Sample.find(aqid)
