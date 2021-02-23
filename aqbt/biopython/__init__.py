"""biopython.py.

Helper functions for BioPython
"""
import hashlib
import itertools
import random
import re
import tempfile
from copy import deepcopy
from itertools import chain
from typing import Iterable
from typing import List
from typing import Tuple
from typing import Union
from uuid import uuid4

import inflection
import networkx as nx
from BCBio import GFF
from Bio import Restriction
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation
from Bio.SeqFeature import ExactPosition
from Bio.SeqFeature import FeatureLocation
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord
from primer3plus.design import primer3

from aqbt.sequence import anneal
from aqbt.sequence import random_sequence
from aqbt.sequence import rc
from aqbt.utils import format_float
from aqbt.utils import random_color
from aqbt.utils import random_slices
from aqbt.utils import sort_cycle
from aqbt.utils.region import Span

# constants
FWD_COLOR = "ApEinfo_fwdcolor"
REV_COLOR = "ApEinfo_revcolor"
CIRCULAR = "circular"
CIRCULAR = "circular"
LINEAR = "linear"
TOPOLOGY = "topology"
DEFAULT_ANNOTATION_TYPE = "misc_feature"


def get_feature_name(feature: SeqFeature):
    return get_feature_qualifiers(feature, "label")


def set_feature_name(feature: SeqFeature, name: str):
    set_feature_qualifier(feature, "label", name)
    return feature


def set_feature_qualifier(feature: SeqFeature, key: str, value):
    feature.qualifiers[key] = [value]


def get_feature_qualifiers(feature: SeqFeature, key: str):
    return feature.qualifiers.get(key, list())[0]


def set_feature_color(feature: SeqFeature, color: str = None, rev_color: str = None):
    if color is None:
        color = random_color()

    if not rev_color:
        rev_color = random_color()

    set_feature_qualifier(feature, FWD_COLOR, color)
    set_feature_qualifier(feature, REV_COLOR, rev_color)


def new_location(i: int, j: int, strand: int) -> FeatureLocation:
    return FeatureLocation(ExactPosition(i), ExactPosition(j), strand=strand)


def new_compound_location(
    indices: List[Union[Tuple[int, int], Tuple[int, int, int]]], strand: int
) -> CompoundLocation:
    locations = []
    for index in indices:
        if not isinstance(index, Tuple):
            raise ValueError(
                "Expects a tuple of integers size 2 or 3, not a {}".format(
                    indices.__class__
                )
            )
        if not len(index) in [2, 3]:
            raise ValueError("Expects a tuple of integers of size 2 or 3")
        if len(index) == 2:
            i, j = index
            s = strand
        elif len(index) == 3:
            i, j, s = index
        else:
            raise ValueError("Must be tuple of 2 or 3 integers")
        if not isinstance(i, int) or not isinstance(j, int) or not isinstance(s, int):
            raise ValueError(
                "Expects a tuple of integers of size 2 or 3. Found {}".format(index)
            )
        locations.append(FeatureLocation(ExactPosition(i), ExactPosition(j), strand=s))
    return CompoundLocation(locations)


def remove_duplicate_features(record: SeqRecord):
    """Remove redundant features."""
    features_by_hash = {}
    for f in record.features:
        fhash = "*".join([str(f.qualifiers["label"]), f.type, str(f.location)])
        features_by_hash[fhash] = f
    record.features = list(features_by_hash.values())
    return record


def new_compound_feature(
    name: str,
    indices: List[Union[Tuple[int, int], Tuple[int, int, int]]],
    strand: int,
    feature_type: str = None,
):
    location = new_compound_location(indices, strand=strand)
    return new_feature(name, location=location, feature_type=feature_type)


def feature_qualifiers_to_snake_case(feature: SeqFeature) -> SeqFeature:
    """Convert the attribute names of feature qualifiers to snake case (lower-
    case, underscore)

    :param feature: BioPython SeqFeature instance
    :return: the same BioPython SeqFeature instance with qualifiers as a new dictionary
    """
    feature.qualifiers = {
        inflection.underscore(k): v for k, v in feature.qualifiers.items()
    }
    return feature


def new_feature(
    name: str,
    location: Union[None, CompoundLocation, FeatureLocation],
    color: str = None,
    feature_type: str = None,
):
    feature = SeqFeature(location=location, qualifiers={}, type=feature_type)
    set_feature_name(feature, name)
    set_feature_color(feature, color=color)
    return feature


def _annotate_feature(
    length: int,
    name: str,
    i: int = None,
    j: int = None,
    cyclic: bool = False,
    feature_type: str = None,
):
    if i is None:
        i = 0
    if j is None:
        j = length

    if cyclic and (j > length or (i > j)):
        if j > length:
            j = j - length
        feature = new_compound_feature(
            name=name,
            indices=[(i, length), (0, j)],
            strand=1,
            feature_type=feature_type,
        )
    else:
        feature = new_feature(
            name=name,
            location=FeatureLocation(ExactPosition(i), ExactPosition(j), strand=1),
            feature_type=feature_type,
        )
    return feature


def annotate(
    record: SeqRecord,
    name: str,
    i: int = None,
    j: int = None,
    cyclic: bool = False,
    annotation_type: str = None,
):
    """Annotate a SeqRecord."""
    if not name:
        raise ValueError("Cannot annotate record with no name.")
    if annotation_type is None:
        annotation_type = DEFAULT_ANNOTATION_TYPE
    feature = _annotate_feature(
        len(record.seq), name, i, j, cyclic, feature_type=annotation_type
    )
    record.features.append(feature)
    return feature


def set_topology(records: List[SeqRecord], topology: str):
    assert topology in [CIRCULAR, LINEAR]
    for r in records:
        r.annotations[TOPOLOGY] = topology


def make_cyclic(records: List[SeqRecord]):
    if isinstance(records, SeqRecord):
        records = [records]
    for r in records:
        r.annotations[TOPOLOGY] = CIRCULAR


def make_linear(records: List[SeqRecord]):
    if isinstance(records, SeqRecord):
        records = [records]
    for r in records:
        r.annotations[TOPOLOGY] = LINEAR


def is_circular(record: SeqRecord):
    return record.annotations.get(TOPOLOGY, False) == CIRCULAR


def is_linear(record: SeqRecord):
    return not is_circular(record)


def new_sequence(
    seqstr: str,
    name: str = None,
    auto_annotate: bool = False,
    cyclic: bool = False,
    annotation_type: str = None,
):
    record = SeqRecord(Seq(seqstr))
    if auto_annotate:
        annotate(record, name, annotation_type=annotation_type)
    if name:
        record.name = name
        record.id = name
        record.description = name

    if cyclic:
        make_cyclic([record])
    else:
        make_linear([record])
    return record


# def slice_feature(feature: SeqFeature, i: int = None, j: int = None):
#     if isinstance(feature)


def _format_slice(slc: slice):
    return "[{}:{}]".format(slc.start, slc.stop)


def slice_with_features(seq: SeqRecord, slc: slice):
    i = slc.start
    j = slc.stop
    if slc.step:
        raise ValueError("Step is not supported.")
    # new_features = []
    # for feature in seq.features:
    #     name = get_feature_name(feature)
    #     new_parts = []
    #     for part in feature.location.parts:
    #         if i in part and j in part:
    #             new_parts.append(new_location(i, j, part.strand))
    #         elif i in part:
    #             new_parts.append(new_location(i, part.end.position, part.strand))
    #         elif j in part:
    #             new_parts.append(new_location(part.start.position, j, part.strand))
    #     if len(new_parts) == 1:
    #         new_feature = deepcopy(feature)
    #         new_feature.location = new_parts[0]
    #         set_feature_name(new_feature, name + _format_slice(slc))
    #         new_features.append(new_feature)
    #     elif len(new_parts) > 2:
    #         new_feature = deepcopy(feature)
    #         new_feature.location = CompoundLocation(new_parts)
    #         set_feature_name(new_feature, name + _format_slice(slc))
    #         new_features.append(new_feature)
    new_seq = seq[i:j]
    # new_seq.features += new_features
    return new_seq


def random_record(
    length: int, name: str = None, auto_annotate: bool = True, cyclic: bool = False
):
    if name is None:
        name = str(uuid4())
    return new_sequence(
        random_sequence(length), name=name, auto_annotate=auto_annotate, cyclic=cyclic
    )


def randomly_annotate(
    rec: SeqRecord,
    feature_length_range: Tuple[int, int],
    feature_name_list: List[str] = None,
) -> SeqRecord:
    for i, j in random_slices(
        feature_length_range[0], feature_length_range[1], len(rec.seq)
    ):
        if not feature_name_list:
            random_name = str(uuid4())[-4:]
        else:
            random_name = random.sample(feature_name_list, k=1)[0]
        annotate(rec, random_name, i, j)
    return rec


def pcr_amplify(
    primers: Iterable[SeqRecord],
    template: SeqRecord,
    cyclic: bool,
    cyclic_buffer: int = 100,
    name: str = None,
    return_matches: bool = False,
    annotate_product: bool = False,
    annotate_primers: bool = False,
) -> Union[List[SeqRecord], Tuple[List[SeqRecord], List[dict], List[dict]]]:
    """

    :param primers: list of primers to use for the amplification. Typically, a
        tuple of two SeqRecords (fwd and rev)
    :param template: template sequence as a SeqRecord
    :param cyclic: whether to treat the template sequence as a cyclic sequence
    :param cyclic_buffer: number of bases to consider for the wrap around origin
        sequence. This buffer allows primers
        to bind to the sequences near the origin. (default: 100)
    :param name: optional name of the new amplicons
    :param return_matches: if True, also return information on the forward and
        reverse matches (default: False)
    :param annotate_product: whether to annotate the final product with the provided name
    :param annotate_primers: whether to annotate the final product with the primer
        names at their binding sites
    :return: either list of amplicons as SeqRecords. If return_matches == True, also
        return the forward and reverse
        primer matches.
    """
    assert isinstance(primers, list) or isinstance(primers, tuple)
    # template_name = template.name or template.id or template.description
    original_template = template
    if cyclic:
        template = template + template + template[:cyclic_buffer]
    fwd_matches, rev_matches = anneal(str(template.seq), [str(x.seq) for x in primers])

    products_by_sequence = {}
    for f, r in itertools.product(fwd_matches, rev_matches):
        i = f["top_strand_slice"][0]
        j = r["top_strand_slice"][1]

        try:
            span = Span(
                i, j, len(original_template.seq), cyclic=cyclic, ignore_wrap=True
            )
        except IndexError as e:
            if not cyclic:
                continue
            else:
                raise e

        fwd = primers[f["primer_index"]]
        rev = primers[r["primer_index"]]
        f_rec = new_sequence(
            f["overhang"], name=fwd.name + "_overhang", auto_annotate=annotate_primers
        )
        r_rec = new_sequence(
            rc(r["overhang"]),
            name=rev.name + "_overhang",
            auto_annotate=annotate_primers,
        )

        product = span.get_slice(original_template)
        # annotate(template_chunk, name="source: {}".format(template_name[:40]))
        if len(f_rec.seq):
            product = f_rec + product
        if len(r_rec.seq):
            product = product + r_rec

        if not name:
            name = original_template.name or original_template.id
            name += "[{}:{}]".format(span.a, span.b)
        if annotate_product:
            annotate(product, name=name)
        if len(product) <= len(original_template.seq):
            products_by_sequence[str(product.seq)] = (product, span.a, span.b)
    product_list = list(products_by_sequence.values())
    product_list.sort(key=lambda x: (x[1], len(product)))
    if return_matches:
        return product_list, fwd_matches, rev_matches
    return product_list


def load_glob(*paths: Tuple[str, ...], format: str = None):
    path_iter = chain(*paths)
    records = []
    for path in path_iter:
        records += SeqIO.parse(path, format=format)
    return records


def load_fasta_glob(path: str) -> List[SeqRecord]:
    return load_glob(path, format="fasta")


def load_genbank_glob(path: str) -> List[SeqRecord]:
    return load_glob(path, format="genbank")


def write_tmp_records(records: List[SeqRecord], format: str) -> List[SeqRecord]:
    fd, tmp_path_handle = tempfile.mkstemp(suffix="." + format)
    SeqIO.write(records, tmp_path_handle, format=format)
    return tmp_path_handle


class GibsonAssembler:
    @staticmethod
    def make_hash(s: str):
        return hashlib.sha1(s.encode("utf-8")).hexdigest()

    @classmethod
    def rnode(cls, record: SeqRecord):
        return cls.make_hash(str(record.name + record.seq))

    @classmethod
    def add_edge(
        cls, g: nx.DiGraph, r1: SeqRecord, r2: SeqRecord, tm: float, match: dict
    ):
        n1 = cls.rnode(r1)
        n2 = cls.rnode(r2)
        if n1 not in g:
            g.add_node(n1, record=r1)
        if n2 not in g:
            g.add_node(n2, record=r2)
        g.add_edge(n1, n2, tm=tm, match=match)

    @classmethod
    def interaction_graph(cls, records: List[SeqRecord]) -> nx.DiGraph:
        g = nx.DiGraph()
        pairs = []
        for r1, r2 in itertools.product(records, records):
            if r1 is not r2:
                pairs.append((r1, r2))
                if r2 is not records[0]:
                    pairs.append((r1, r2.reverse_complement(name=r2.name + "_rc")))

        for r1, r2 in pairs:
            fwd, rev = anneal(str(r1.seq), [str(r2.seq)])
            for f in fwd:
                if f["top_strand_slice"][0] == 0:
                    anneal_seq = f["anneal"]
                    tm = primer3.calcTm(anneal_seq[-60:])
                    cls.add_edge(g, r2, r1, tm=tm, match=f)
        return g

    @classmethod
    def _collect_features(cls, records: List[SeqRecord]):
        seq_to_features = {}
        for record in records:
            for feature in record.features:
                s = str(feature.extract(record.seq)).upper()
                if s not in seq_to_features:
                    seq_to_features[s] = feature
        return seq_to_features

    @classmethod
    def _restore_features(cls, seq_to_features, record, cyclic: bool):
        for seq in seq_to_features:
            record_seq = str(record.seq)
            if cyclic:
                record_seq += record_seq
            matches = list(re.finditer(seq, record_seq, re.IGNORECASE))
            if len(matches) >= 1:
                feature = seq_to_features[seq]
                new_feature = deepcopy(feature)
                for match in matches:
                    span = match.span()

                    if span[0] < len(record.seq):
                        if span[1] >= len(record.seq):
                            f = _annotate_feature(
                                len(record.seq), "", span[0], span[1], cyclic=True
                            )
                            new_feature.location = f.location
                        else:
                            new_feature.location = new_location(
                                span[0], span[1], feature.location.strand
                            )
                    record.features.append(new_feature)

    @classmethod
    def assemble_record_from_cycle(
        cls,
        cycle: List[str],
        graph: nx.DiGraph,
        annotate_sources: bool = False,
        annotate_junctions: bool = False,
    ) -> SeqRecord:

        source_indices = []
        i = 0

        records = [graph.nodes[n]["record"] for n in cycle]
        stored_features = cls._collect_features(records)

        record = SeqRecord(Seq(""))

        c1 = cycle[-1:] + cycle[:-1]
        c2 = cycle
        c3 = cycle[1:] + cycle[:1]

        for n1, n2, n3 in zip(c1, c2, c3):
            r1, r2, r3 = [graph.nodes[n]["record"] for n in [n1, n2, n3]]

            edata12 = graph.get_edge_data(n1, n2)
            edata23 = graph.get_edge_data(n2, n3)

            s12 = edata12["match"]["top_strand_slice"]
            s23 = edata23["match"]["top_strand_slice"]

            o12 = r2[s12[0] : s12[1]]
            o23 = r3[s23[0] : s23[1]]

            trunc_r2 = slice_with_features(r2, slice(len(o12), -len(o23)))

            if annotate_junctions:
                annotate(o12, "junction: tm {}".format(format_float(edata12["tm"])))
            record += o12
            record += trunc_r2

            source_indices.append((i, i + len(r2)))
            i = i + len(r2) - len(o23)

        for i, (a, b) in enumerate(source_indices):
            if annotate_sources:
                annotate(
                    record,
                    "source: fragment {} ({})".format(i, r2.name),
                    a,
                    b,
                    cyclic=True,
                )

        cls._restore_features(stored_features, record, cyclic=True)
        remove_duplicate_features(record)
        return record

    @classmethod
    def make_cyclic_assemblies(
        cls,
        records: List[SeqRecord],
        annotate_sources: bool = False,
        annotate_junctions: bool = False,
    ) -> List[SeqRecord]:
        g = cls.interaction_graph(records)

        node_rank = {}
        for i, r in enumerate(records):
            node_rank[cls.rnode(r)] = i

        cycles = list(nx.simple_cycles(g))
        cycles = [sort_cycle(cycle, key=lambda n: node_rank[n]) for cycle in cycles]

        assembled_records = []
        for cycle in cycles:
            record = cls.assemble_record_from_cycle(
                cycle,
                g,
                annotate_sources=annotate_sources,
                annotate_junctions=annotate_junctions,
            )
            assembled_records.append(record)
        make_cyclic(assembled_records)
        return assembled_records


def digest(record: SeqRecord, restriction_enzyme: str):
    """Digest a SeqRecord with an enzyme.

    :param record: record to digest
    :param restriction_enzyme: name of restriction enzyme
    :return: list of records
    """
    enzyme = getattr(Restriction, restriction_enzyme)
    linear = is_linear(record)
    indices = enzyme.search(record.seq, linear=linear)
    indices = [i - 1 for i in indices]
    pairs = zip(indices[:-1], indices[1:])
    records = [record]
    for i1, i2 in pairs:
        records.append(record[i1:i2])
    if not linear:
        records.append(record[indices[-1] :] + record[: indices[0]])
    return records


def _feature_qualifiers_to_snake_case(feature: SeqFeature) -> SeqFeature:
    feature.qualifiers = {
        inflection.underscore(k): v for k, v in feature.qualifiers.items()
    }
    return feature


def _correct_phase(feature: SeqFeature):
    """Corrects the 'phase' attribute for SeqFeatures. Some parsers produce
    invalid phases (must be 0, 1, or 2). Invalid phases are delected.

    :param feature:
    :return:
    """
    phase_key = "Phase"
    if phase_key in feature.qualifiers:
        if feature.qualifiers[phase_key] not in (0, 1, 2):
            del feature.qualifiers[phase_key]


def prepare_record_for_gff(record: SeqRecord, in_place: bool = True):
    """Prepare a SeqRecord for GFF conversion.

    1. Converts feature.qualifier['label'] to 'ID' and 'Name'
    :param record:
    :param in_place:
    :return:
    """
    if not in_place:
        record = deepcopy(record)

    for feature in record.features:
        _correct_phase(feature)

    for feature in record.features:

        label = feature.qualifiers.get("label", [None])[0]
        if label:
            if "ID" not in feature.qualifiers:
                feature.qualifiers["ID"] = label
            if "Name" not in feature.qualifiers:
                feature.qualifiers["Name"] = feature.qualifiers["ID"]

        _feature_qualifiers_to_snake_case(feature)

    return record


def to_gff(
    records: List[SeqRecord],
    out_path: str,
    include_fasta: bool = True,
    in_place: bool = False,
):
    """Convert a list of SeqRecords to an annotated GFF file.

    :param records: list of SeqRecords.
    :param out_path: where to save the GFF file.
    :param include_fasta: include the FASTA sequence in the GFF file
    :param in_place: Whether to alter the list of records.
    :return: out_path
    """
    to_gff_recs = []
    for rec in records:
        to_gff_recs.append(prepare_record_for_gff(rec, in_place=in_place))

    with open(out_path, "w") as f:
        GFF.write(to_gff_recs, f, include_fasta=include_fasta)
    return out_path


def from_gffs(paths: List[str]) -> List[SeqRecord]:
    records = list(GFF.parse(paths))
    return records


make_cyclic_assemblies = GibsonAssembler.make_cyclic_assemblies
