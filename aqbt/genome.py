import itertools
from copy import deepcopy
from typing import Any
from typing import Callable
from typing import Dict
from typing import List
from typing import Tuple

import networkx as nx
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from dasi.utils.region import Span
from pyblast import BioBlast
from pydent import AqSession
from pydent.models import Sample

from aqbt import biopython
from aqbt.logger import logger


class IntegrationError(Exception):
    """Generic exceptions for integration events."""


def _generate_integration_graph(
    sess: AqSession,
    sample: Sample,
    field_names: Tuple[str, ...] = ("Integrant", "Parent"),
) -> nx.DiGraph:
    """Generate a DAG from the integrant and parent of the provided aquarium
    sample.

    :param sess: AqSession instance
    :param sample: Aquarium sample
    :param field_names: list of field names to look for to generate the DAG. default: (Integrant, Parent)
    :return:
    """
    g = sess.browser.sample_network([sample])
    new_graph = nx.DiGraph()
    for n1, n2, edata in g.edges(data=True):
        if edata["attr_dict"]["field_value_name"] in field_names:
            new_graph.add_edge(n1, n2, **edata)
    return new_graph


# def build_edge_genome_from_aq_sample(
#     edge_api: EdgeAPI,
#     registry: KlavinsLabRegistry,
#     sample: Sample,
#     aq_id_to_edge_id: Dict[int, int],
#     restriction_enzyme: str,
#     genome_prefix: str = "",
#     genome_suffix: str = "",
# ):
#     """
#
#     :param edge_api: edge API instance
#     :param registry: KlavinsLabRegistry instance
#     :param sample: aquarium sample to produce GFF file
#     :param aq_id_to_edge_id: base strain aq_id to edge_id (e.g. convert Aquarium 22800 to edge id 1 {22800: 1})
#     :param restriction_enzyme: name of the restriction enzyme to linearize the plasmid.
#     :return:
#     """
#     session = registry.session
#     if restriction_enzyme:
#         enzyme = getattr(Restriction, restriction_enzyme)
#     else:
#         enzyme = None
#     with session.with_cache(timeout=60) as sess:
#         new_graph = _generate_integration_graph(sess, sample)
#         yeast_type = sess.SampleType.find_by_name("Yeast Strain")
#         nodes = nx.topological_sort(new_graph)
#
#         for ntype, nid in nodes:
#             edge_id = aq_id_to_edge_id.get(nid, None)
#             yeast = sess.Sample.find(nid)
#             if yeast.sample_type_id == yeast_type.id:
#                 if not edge_id:
#                     plasmid = yeast.properties["Integrant"]
#                     parent = yeast.properties["Parent"]
#
#                     plasmid_seq = registry.get_sequence(plasmid)
#                     print(plasmid_seq.web_url)
#
#                     plasmid_record = convert(plasmid_seq, to="SeqRecord")
#                     plasmid_is_linear = not is_circular(plasmid_record)
#
#                     if enzyme:
#
#                         digestion_parts = enzyme.catalyse(
#                             plasmid_record.seq, linear=plasmid_is_linear
#                         )
#                     else:
#                         digestion_parts = [plasmid_record.seq]
#
#                     ok_parts = []
#                     genome_name = (
#                         genome_prefix
#                         + str(yeast.id)
#                         + "__"
#                         + yeast.name
#                         + genome_suffix
#                     )
#                     recomb = partial(
#                         edge_api.recombine,
#                         parent.edge_id,
#                         genome_name=genome_name,
#                         homology_arm_length=200,
#                         cassette_name="UWBF{}_{}".format(plasmid.id, plasmid.name),
#                         notes="Linearized with {RE}".format(RE=restriction_enzyme),
#                     )
#
#                     for part in digestion_parts:
#                         response = recomb(
#                             create=False, cassette=str(part), verbose=True
#                         )
#                         print(response)
#                         if isinstance(response, list) and response:
#                             ok_parts.append(part)
#
#                     if len(ok_parts) > 1:
#                         for part in ok_parts:
#                             print(part)
#                         raise ValueError(
#                             "More than one part integrates for {}".format(plasmid.name)
#                         )
#                     elif len(ok_parts) == 0:
#                         raise ValueError(
#                             "No parts integrate for {}".format(plasmid.name)
#                         )
#
#                     response = recomb(create=True, cassette=str(ok_parts[0]))
#                     assert response["id"]
#                     yeast.edge_id = response["id"]
#                 else:
#                     yeast.edge_id = edge_id

# TODO: set appropriate features
def _get_integration_fragments(plasmid_record, restriction_enzyme: str) -> List[Seq]:
    # plasmid_is_linear = not is_circular(plasmid_record)

    # perform (optional) digest
    if restriction_enzyme:
        integration_records = biopython.digest(plasmid_record, restriction_enzyme)
        for i, record in enumerate(integration_records):
            suffix = "_{}_digest_{}".format(restriction_enzyme, i)
            record.id = plasmid_record.id + suffix
            record.name = plasmid_record.name + suffix
            record.annotations["origin_id"] = plasmid_record.id
    else:
        integration_records = [deepcopy(plasmid_record)]

    for r in integration_records:
        biopython.annotate(r, name="Integration")

    return integration_records


def get_span_info(result_entry: Dict) -> Tuple[int, int, int]:
    """Return the start, stop, and strand information. Start and stop is in
    reference to the top strand.

    e.g.::

        1, 100, 1 means the span goes from 1 to 100 (inclusively) on the top strand

        1, 100, -1 means the span goes from 100 to 1 (inclusively) on the bottom strand.

    :param result_entry:
    :return:
    """
    a = result_entry["start"]
    b = result_entry["end"]
    strand = result_entry["strand"]
    if strand == -1:
        a, b = a, b
    return a, b, strand


def get_strand(result: Dict) -> int:
    strand = result["query"]["strand"] * result["subject"]["strand"]
    if strand not in [1, -1]:
        raise ValueError(r"Strand must be in (1, -1), not {}".format(strand))
    return strand


def simplify_result(result: Dict) -> Tuple[int, int, int]:
    strand = get_strand(result)
    query_span = get_span_info(result["query"])
    subject_span = get_span_info()
    assert subject_span[0] < subject_span[1]
    assert query_span[0] < query_span[1]
    return {
        "query": (query_span[0], query_span[1]),
        "subject": (subject_span[0], subject_span[1]),
        "strand": strand,
    }


# class Alignment(object):
#
#     def __init__(self, locus_start, locus_end, cassette_start, cassette_end, direction):
#         """
#
#         :param locus_start: inslusive start of homology (top strand)
#         :param locus_end: exclusive end of homology (top strand)
#         :param cassette_start: inclusive start of cassette (top strand)
#         :param cassette_end: exclusive end of cassette (top strand)
#         :param direction: 1 or -1.
#             "1" means cassette and locus are pointing in same direction.
#             "-1" means cassette and locus are on different strands.
#         """
#         assert locus_start <= locus_end
#         assert cassette_start <= cassette_end
#         self.locus_start = locus_start
#         self.locus_end = locus_end
#         self.cassette_start = cassette_start
#         self.cassette_end = cassette_end
#         self.direction = direction
#
#
# def _correct_indices(blast_start, blast_end, blast_strand):
#     if blast_strand == 1:
#         a, b, s = blast_start, blast_end, blast_strand
#     elif blast_strand == -1:
#         a, b, s = blast_end, blast_start, blast_strand
#     # inclusive to exclusive
#     b += 1
#     return a, b, s
#
# def blast_res_to_alignment(result):
#
#     def f(e):
#         return _correct_indices(e['start'], e['end'], e['strand'])
#     s = f(result['subject'])
#     q = f(result['query'])
#     direction = s[2] * q[2]
#     return Alignment(s[0], s[1], q[0], q[1], direction)
#
# def is_valid_integration(a1: Alignment, a2: Alignment, max_deletion_length):
#     error_msgs = []
#     if a1.direction != a2.direction:
#         error_msgs.append("Homology regions are on opposite strands.")
#
#     query_direction = a1.cassette_end <= a2.cassette_start
#     subject_direction = a1.locus_start <= a2.locus_end
#     if not query_direction == subject_direction:
#         error_msgs.append("Homology regions are swapped.")
#
#     region_to_delete = a2.locus_start - a1.locus_end
#     cassette_region = a2.cassette_start - a1.cassette_end
#
#     return {
#         'is_valid': len(error_msgs) == 0,
#         'errors': error_msgs
#     }


class Alignment:
    def __init__(self, subject: Span, query: Span, direction: int):
        self.subject = subject
        self.query = query
        self.direction = direction


def group_by(arr, group_by_func: Callable) -> Dict[Any, List]:
    grouped = {}
    for x in arr:
        key = group_by_func(x)
        grouped.setdefault(key, list())
        grouped[key].append(x)
    return grouped


def get_left_flank_span(span: Span) -> Span:
    return span.new(0, span.a)


def get_right_flank_span(span: Span) -> Span:
    if span.b == span.context_length:
        return span.new(0, 0)
    return span.new(span.b, span.context_length)


def get_flanking_spans(span: Span) -> Tuple[Span, Span]:
    return get_left_flank_span(span), get_right_flank_span(span)


def parse_result_to_span(data, inclusive=True, input_index=1, output_index=None):
    if output_index is None:
        output_index = input_index

    s, e, length = data["start"], data["raw_end"], data["origin_sequence_length"]
    if data["strand"] == -1:
        s, e = e, s
    c = data["circular"]
    if inclusive:
        e += 1
    span = Span(s, e, length, cyclic=c, ignore_wrap=False, index=input_index)
    if input_index != output_index:
        span = span.reindex(output_index)
    return span


class IntegrationSite:
    def __init__(
        self,
        chr: SeqRecord,
        integrative_dna: SeqRecord,
        direction: int,
        chr_five_prime_flanking_region_span: Span,
        chr_three_prime_flanking_region_span: Span,
        chr_excision_region_span: Span,
        chr_five_prime_homology_arm_span: Span,
        chr_three_prime_homology_arm_span: Span,
        chr_integration_cassette_span: Span,
        insert_five_prime_homology_arm_span: Span,
        insert_three_prime_homology_arm_span: Span,
        insert_integration_cassette_span: Span,
    ):
        self.chr = chr
        self.integrative_dna = integrative_dna
        self.direction = direction
        self.chr_five_prime_flanking_region_span = chr_five_prime_flanking_region_span
        self.chr_three_prime_flanking_region_span = chr_three_prime_flanking_region_span
        self.chr_excision_region_span = chr_excision_region_span
        self.chr_five_prime_homology_arm_span = chr_five_prime_homology_arm_span
        self.chr_three_prime_homology_arm_span = chr_three_prime_homology_arm_span
        self.chr_integration_cassette_span = chr_integration_cassette_span
        self.insert_five_prime_homology_arm_span = insert_five_prime_homology_arm_span
        self.insert_three_prime_homology_arm_span = insert_three_prime_homology_arm_span
        self.insert_integration_cassette_span = insert_integration_cassette_span


class GenomeIntegrator:
    def __init__(
        self,
        e_value_threshold: float = 1e-10,
        minimum_homology_length: int = 30,
        max_genomic_deletion: int = 5000,
        max_integration_overhang: int = 30,
    ):
        self.e_value_threshold = e_value_threshold
        self.min_homology_length = minimum_homology_length
        self.max_genomic_deletion = max_genomic_deletion
        self.max_integration_overhang = max_integration_overhang
        self.logger = logger(self)

    def is_valid_integration(self, left: Alignment, right: Alignment):
        """

        :param left: Left most alignment (in relation to integration cassette)
        :param right: Right most alignment (in relation to integration cassette)
        :return:
        """

        self.logger.debug(
            "Integration Validation:\n"
            "Left subject:  {}\n"
            "Right subject: {}\n"
            "Left query:   {}\n"
            "Right query:  {}\n"
            "Direction:    {}".format(
                left.subject, right.subject, left.query, right.query, left.direction
            )
        )

        errors = []
        if not left.direction == right.direction:
            errors.append(
                "left and right alignments are pointing in different directions"
            )
            self.logger.debug(errors[-1])

        # if a1.query.b < a2.query.a:
        #     left, right = a1, a2
        # elif a2.query.b < a1.query.a:
        #     left, right = a2, a1

        if left.direction == -1:
            a = right.subject.a
            b = left.subject.b
        else:
            a = left.subject.b
            b = right.subject.a
        if a == b:
            print()

        try:
            gap = right.subject.new(a, b)
        except IndexError:
            gap = right.subject.intersection(left.subject)

        if gap is None:
            errors.append("No gap or intersection between homologies found.")
            self.logger.debug(errors[-1])
            return {"errors": errors}

        left_overhang_span = get_left_flank_span(left.query)
        right_overhang_span = get_right_flank_span(right.query)

        if len(left_overhang_span) > self.max_integration_overhang:
            errors.append(
                "left_overhang_span {} > {}".format(
                    len(left_overhang_span), self.max_integration_overhang
                )
            )
            self.logger.debug(errors[-1])

        if len(right_overhang_span) > self.max_integration_overhang:
            errors.append(
                "right_overhang_span {} > {}".format(
                    len(right_overhang_span), self.max_integration_overhang
                )
            )
            self.logger.debug(errors[-1])

        if len(gap) > self.max_genomic_deletion:
            errors.append(
                "len(gap) {} > {}".format(len(gap), self.max_genomic_deletion)
            )
            self.logger.debug(errors[-1])
        if errors:
            return {"errors": errors}
        self.logger.debug("Integration passed validation!")
        return {"left": left, "right": right, "gap": gap}

    def _group_possible_integration_loci(
        self, blast_results: Dict, group_by_func: Callable = None
    ) -> Dict[Tuple[str, str, int], List[Dict[str, Any]]]:
        """Given a list of subject SeqRecords and a list of query SeqRecords,
        group *possible* regions an integration cassette could be.

        This groups blast results by their origin id and direction.

        :param subjects: list of subjects (e.g. list of chromosomes)
        :param queries: list of queries (e.g. list of integration cassettes)
        :param e_value_threshold: the largest e value to consider a result
        :param min_homology_length: minimum alignment length to consider
        :param group_by_func: optional group_by_function for each blast result.
        :return: grouped blast results
        """

        valid_results = []
        for result in blast_results:
            if (
                result["meta"]["evalue"] < self.e_value_threshold
                and result["meta"]["alignment length"] >= self.min_homology_length
            ):
                valid_results.append(result)

        if group_by_func is None:

            def group_by_func(r):
                key1 = r["query"]["origin_key"]
                key2 = r["subject"]["origin_key"]
                key3 = r["subject"]["strand"]
                key = (key1, key2, key3)
                return key

        return group_by(blast_results, group_by_func)

    def get_valid_integration_sites(
        self, subjects: List[SeqRecord], integrants: List[SeqRecord]
    ) -> List[IntegrationSite]:
        """Generate a list of valid integration sites.

        Result keys:

        .. code-block:: python

            {
                # the region from the beginning of the chromosome to
                # the 5' homology region
                "chr_five_prime_flanking_region_span"

                # the region from the end of the 3' homology region to
                # the end of the chromosome
                "chr_three_prime_flanking_region_span"

                # the region of the original chromosome that was removed
                # by the integration
                "chr_excision_region_span"

                # the region of the 5' homology
                "chr_five_prime_homology_arm_span"

                # the region of the 3' homology
                "chr_three_prime_homology_arm_span"

                # the region from the end of the 5' homology to the start
                # of the 3' homology
                "integration_cassette_span"

                # the region on the integration cassette, including the homology
                # arms, on the integrative plasmid or fragment
                "cassette_span"

                # SeqRecord of the cassette
                "cassette"

                # id of the chromosome
                "chr_id"

                # direction of the integration
                "direction"

                # SeqRecord of the left flank
                "left_flank"

                # SeqRecord of the right flank
                "right_flank"
            }

        :param subjects:
        :param integrants:
        :return:
        """
        blast = BioBlast(subjects, integrants)
        blast_results = blast.blastn()

        grouped = self._group_possible_integration_loci(blast_results)

        def entry_to_span(entry) -> Span:
            return parse_result_to_span(entry, input_index=1, output_index=0)

        def blast_res_to_alignment(result):
            return Alignment(
                subject=entry_to_span(result["subject"]),
                query=entry_to_span(result["query"]),
                direction=result["subject"]["strand"],
            )

        integration_sites = []

        for k, regions in grouped.items():
            if len(regions) > 1:
                for left_result, right_result in itertools.combinations(regions, r=2):
                    left = blast_res_to_alignment(left_result)
                    right = blast_res_to_alignment(right_result)
                    if left.query.b > right.query.a:
                        left, right = right, left
                        left_result, right_result = right_result, left_result

                    is_valid = self.is_valid_integration(left=left, right=right)
                    if is_valid and "errors" not in is_valid:
                        if left.direction == 1:
                            left_flank_span = get_left_flank_span(left.subject)
                            right_flank_span = get_right_flank_span(right.subject)
                            left_hom_span = left.subject
                            right_hom_span = right.subject
                        elif left.direction == -1:
                            left_flank_span = get_left_flank_span(right.subject)
                            right_flank_span = get_right_flank_span(left.subject)
                            left_hom_span = right.subject
                            right_hom_span = left.subject

                        else:
                            raise ValueError("Direction not understood.")
                        gap_span = is_valid["gap"]
                        cassette_span = left.query.new(left.query.a, right.query.b)
                        chr = blast.seq_db.get_origin(
                            left_result["subject"]["origin_key"]
                        )
                        cassette = blast.seq_db.get_origin(
                            left_result["query"]["origin_key"]
                        )
                        cassette.name = cassette.name + "[{}:{}]".format(
                            cassette_span.a, cassette_span.b
                        )
                        # left_flank_span.get_slice = Span.get_slice
                        # right_flank_span.get_slice = Span.get_slice
                        # cassette_span.get_slice = Span.get_slice

                        integration_sites.append(
                            IntegrationSite(
                                chr=chr,
                                integrative_dna=cassette,
                                direction=left.direction,
                                chr_five_prime_flanking_region_span=left_flank_span,
                                chr_three_prime_flanking_region_span=right_flank_span,
                                chr_excision_region_span=gap_span,
                                chr_five_prime_homology_arm_span=left_hom_span,
                                chr_three_prime_homology_arm_span=right_hom_span,
                                chr_integration_cassette_span=left.subject.new(
                                    left_flank_span.b, right_flank_span.a
                                ),
                                insert_five_prime_homology_arm_span=left.query,
                                insert_three_prime_homology_arm_span=right.query,
                                insert_integration_cassette_span=left.query.new(
                                    left.query.b, right.query.a
                                ),
                            )
                        )
                        # integration_sites.append(
                        #     {
                        #         # the region from the beginning of the chromosome to
                        #         # the 5' homology region
                        #         "chr_five_prime_flanking_region_span": left_flank_span,
                        #
                        #         # the region from the end of the 3' homology region to
                        #         # the end of the chromosome
                        #         "chr_three_prime_flanking_region_span": right_flank_span,
                        #
                        #         # the region of the original chromosome that was removed
                        #         # by the integration
                        #         "chr_excision_region_span": gap_span,
                        #
                        #         # the region of the 5' homology
                        #         "chr_five_prime_homology_arm_span": left.subject,
                        #
                        #         # the region of the 3' homology
                        #         "chr_three_prime_homology_arm_span": right.subject,
                        #
                        #         # the region from the end of the 5' homology to the start
                        #         # of the 3' homology
                        #         "chr_integration_cassette_span": left.subject.new(left.query.b, right.query.a),
                        #
                        #         # the region on the integrative DNA from the end of its
                        #         # 5' homology arm to the start of the 3' homology region
                        #         "insert_integration_cassette_span": left.query.new(left.query.b, right.query.a),
                        #
                        #         # the region on the integrative DNA from the beginning of its
                        #         # 5' homology region to the end of the 5' homology region
                        #         "insert_five_prime_homology_arm_span": left.query.copy(),
                        #
                        #         # the region on the integrative DNA from the beginning of its
                        #         # 3' homology region to the end of the 3' homology region
                        #         "insert_three_prime_homology_arm_span": right.query.copy(),
                        #
                        #         # the region on the integration cassette, including the homology
                        #         # arms, on the integrative plasmid or fragment
                        #         "cassette_span": cassette_span,
                        #
                        #         # SeqRecord of the cassette
                        #         "cassette": cassette_span.get_slice(cassette),
                        #
                        #         # id of the chromosome
                        #         "chr_id": blast.seq_db.get_origin(
                        #             left_result["subject"]["origin_key"]
                        #         ).id,
                        #
                        #         # direction of the integration
                        #         "direction": left.direction,
                        #
                        #         # SeqRecord of the left flank
                        #         "left_flank": left_flank_span.get_slice(chr),
                        #
                        #         # SeqRecord of the right flank
                        #         "right_flank": right_flank_span.get_slice(chr),
                        #     }
                        # )

        return integration_sites

    def integrate(
        self,
        chromosomes: List[SeqRecord],
        integrants: List[SeqRecord],
    ) -> List[Dict[str, Any]]:
        """From a list of chromosomes and list of integrants, generate a new
        engineered genome.

        :param chromosomes: list of chromosomes as SeqRecords
        :param integrants: list of integrants as SeqRecords
        :return: engineered genome as a dictionary
        """
        sites = self.get_valid_integration_sites(chromosomes, integrants)

        result = {"genome": None, "integration": None, "errors": [], "note": None}

        if len(sites) > 1:
            err_msg = "More than one integration site found."
            result["errors"].append(err_msg)
            self.logger.error(err_msg)
        if len(sites) == 0:
            result["errors"].append("No integration sites found.")
            return result
        site = sites[0]

        new_genome = []
        for chr in chromosomes:
            if chr.id == site.chr.id:
                cassette = site.integrative_dna
                if site.direction == -1:
                    cassette = cassette.reverse_complement()
                left_flank = site.chr_five_prime_flanking_region_span.get_slice(chr)
                left_hom = site.insert_five_prime_homology_arm_span.get_slice(cassette)
                insert = site.insert_integration_cassette_span.get_slice(cassette)
                right_hom = site.insert_three_prime_homology_arm_span.get_slice(
                    cassette
                )
                right_flank = site.chr_three_prime_flanking_region_span.get_slice(chr)

                new_chr = left_flank + left_hom + insert + right_hom + right_flank

                new_chr.id = chr.id
                new_chr.name = chr.name
                new_genome.append(new_chr)
            else:
                new_genome.append(chr)

        result["genome"] = new_genome
        result["integration"] = site
        result[
            "note"
        ] = "{len}bp integration cassette integrated into {chr} at {locus} (direction={dir}).'".format(
            len=len(site.integrative_dna),
            chr=site.chr.id,
            locus=(site.chr_excision_region_span.a, site.chr_excision_region_span.b),
            dir=site.direction,
        )
        return result
