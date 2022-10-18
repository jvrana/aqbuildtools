import networkx as nx
from benchlingapi.models import DNASequence
from pydent.models import Sample

from aqbt import biopython
from aqbt.aquarium import pydent_utils
from aqbt.aquarium.pydent_utils import Constants as C
import logging
import typing as t
from Bio.SeqRecord import SeqRecord
logger = logging.getLogger(__file__)

class Resolver:
    def __init__(self, session, registry):
        self.registry = registry.copy()
        self.registry.using_cache = True
        self.plasmid_type = session.SampleType.find_by_name(C.PLASMID)
        self.fragment_type = session.SampleType.find_by_name(C.FRAGMENT)
        self.dna_types = [self.plasmid_type, self.fragment_type]
        self.session = session
        self.resolved_sequences = {}
        self.force_build_at_depths = []
        self.force_build_types = []
        self.force_build_sample_ids = []

    def info(self, msg):
        print(msg)
        logger.info(msg)

    @staticmethod
    def keystone_ops(session, sample):
        ops = pydent_utils.ops_using_samples(session, [sample], role="output")
        return [op for op in ops if pydent_utils.is_keystone_op(op)]

    @staticmethod
    def select_ops(ops):
        if not ops:
            return None
        return ops[0]

    def get_input_samples(self, op):
        samples = []
        for fv in op.inputs:
            if fv.sample:
                if fv.sample.sample_type_id in [
                    self.plasmid_type.id,
                    self.fragment_type.id,
                ]:
                    samples.append(fv.sample)
        return samples

    def get_input_sequences(self, sample, current_depth):
        ops = self.keystone_ops(self.session, sample)
        op = self.select_ops(ops)
        if not op:
            self.info("ERROR: no operations produce '{}'".format(sample.name))
            return {}
        self.session.browser.get([op], {"field_values": {"sample": "field_values"}})
        input_samples = self.get_input_samples(op)
        input_sequences = {}

        for input_sample in input_samples:
            input_sequence = self._resolve(
                sample, input_sample, depth=current_depth + 1
            )
            input_sequences[input_sample.id] = input_sequence
        return input_sequences

    def _try_find(self, sample, current_depth):
        if sample.id in self.resolved_sequences:
            self.info("sequence already resolved")
            return self.resolved_sequences[sample.id]
        if current_depth in self.force_build_at_depths:
            self.info(
                "registry ignored at depth. Depth {} is in {}".format(
                    current_depth, self.force_build_at_depths
                )
            )
            return None
        elif sample.sample_type.name in self.force_build_types:
            self.info(
                "registry ignored. Type {} is in {}".format(
                    sample.sample_type.name, self.force_build_types
                )
            )
            return None
        elif sample.id in self.force_build_sample_ids:
            self.info(
                "registry ignored. Sample id {} is in {}".format(
                    sample.id, self.force_build_sample_ids
                )
            )
            return None
        return self.registry.get_sequence(sample)

    def build_plasmid(self, sample, current_depth):
        self.info('WARNING: No sequence. Looking for how to build "{}"'.format(sample.name))
        input_sequences = self.get_input_sequences(sample, current_depth)

        failed = False
        for sid, input_seq in input_sequences.items():
            if input_seq is None:
                self.info('ERROR: no sequence for "{}"'.format(sid))
                failed = True
        if failed:
            return

        if not input_sequences:
            self.info('ERROR: no input sequences found for "{}"'.format(sample.id))
            return

        input_records = self.registry.connector.convert(
            list(input_sequences.values()), to="SeqRecord"
        )
        assemblies = biopython.make_cyclic_assemblies(input_records)
        if not assemblies:
            self.info("ERROR: no cyclic assembly for '{}'".format(sample.name))
        elif len(assemblies) == 1:
            self.info("SUCCESS: Resolved plasmid from gibson assembly!")
            return self.registry.connector.convert(assemblies[0], to="DNASequence")
        else:
            self.info("ERROR: more than one sequence")

    def resolve_plasmid(self, sample, current_depth) -> t.Union[None, DNASequence]:
        seq = self._try_find(sample, current_depth)
        if not seq:
            seq = self.build_plasmid(sample, current_depth)
        return seq

    def build_fragment(self, sample: Sample) -> t.Union[None, DNASequence]:
        self.info(f"building fragment sample={sample.id}:{sample.name}")

        products = self.make_pcr_product(sample)

        size_select_perc_diff_threshold = 10  # select products within 10% of expected length
        if len(products) > 1:

            # then select the one closest to the given length
            self.info("Found more than one product. Lengths=" + ','.join(str(len(p[0].seq)) for p in products))
            expected_length = sample.properties['Length']
            if not expected_length:
                self.info("No expected length found. Cannot size select multiple products without a length")
                return None
            else:
                size_selected_products = []
                for p in products:
                    diff = len(p[0].seq) - expected_length
                    perc_diff = (diff / expected_length) * 100
                    x = size_select_perc_diff_threshold
                    if perc_diff <= x and perc_diff >= -x:
                        size_selected_products.append(p)
                if len(size_selected_products) > 1:
                    self.info(f"Multiple products found within 10% of expected length {expected_length}bp")
                elif len(size_selected_products) == 0:
                    self.info(f"fNo products found within 10% of expected lenght {expected_length}bp")
                else:
                    dna_product = self.registry.connector.convert(
                        p[0], to="DNASequence"
                    )
                    dna_product.is_circular = False
                    return dna_product
            self.info("ERROR: more than one product")
            return None
        elif len(products) == 1:
            dna_product = self.registry.connector.convert(
                products[0][0], to="DNASequence"
            )
            dna_product.is_circular = False
            self.info("SUCCESS: Build fragment from primers and template.")
            return dna_product
        else:
            self.info("ERROR: no products")
            return None

    def resolve_fragment(self, sample, current_depth: int) -> t.Union[None, DNASequence]:
        seq = self._try_find(sample, current_depth)
        if not seq:
            seq = self.build_fragment(sample)
        return seq

    def make_pcr_product(self, fragment) -> t.List[t.Tuple[SeqRecord, int, int]]:
        fwd_primer_sample = fragment.properties["Forward Primer"]
        rev_primer_sample = fragment.properties["Reverse Primer"]
        template_sample = fragment.properties["Template"]

        if not fwd_primer_sample:
            self.info("ERROR: no 'Forward Primer' found in sample")
            return []
        if not rev_primer_sample:
            self.info("ERROR: no 'Reverse Primer' found in sample")
            return []
        if not template_sample:
            self.info("ERROR: no 'Template' found in sample")
            return []

        fwd = self.registry.get_primer_sequence(fwd_primer_sample)
        rev = self.registry.get_primer_sequence(rev_primer_sample)
        template = self.registry.get_sequence(template_sample)

        if not fwd:
            self.info("ERROR: no fwd primer")
            return []
        if not rev:
            self.info("ERROR: no rev primer")
            return []
        if not template:
            self.info("ERROR: no template")
            return []

        template_record = self.registry.connector.convert(template, to="SeqRecord")

        product_records = biopython.pcr_amplify(
            (fwd, rev), template_record, cyclic=template.is_circular, name=fragment.name
        )
        return product_records

    def _add_edge(self, src, dest):
        pass

    def _resolve(self, source_sample, dest_sample, depth: int) -> t.Union[None, DNASequence]:
        self.info("Resolving sequence for '{}'".format(dest_sample.name))
        self._add_edge(source_sample, dest_sample)
        if dest_sample.sample_type_id == self.plasmid_type.id:
            seq = self.resolve_plasmid(dest_sample, current_depth=depth)
        elif dest_sample.sample_type_id == self.fragment_type.id:
            seq = self.resolve_fragment(dest_sample, current_depth=depth)
        self.resolved_sequences[dest_sample.id] = seq
        if not seq:
            return
        else:
            seq.name = dest_sample.name
            seq.description = dest_sample.description
            if seq.id and "seq_" not in seq.id:
                seq.id = None
            if hasattr(seq, "primers") and seq.primers:
                seq.primers = []
            return seq

    def resolve_sequence(self, sample: Sample) -> t.Union[None, DNASequence]:
        return self._resolve(None, sample, 0)

    def register_new_sequences(self):
        for sid, seq in self.resolved_sequences.items():
            sample = self.session.Sample.find(sid)
            if seq and not self.registry.find_in_cache(sample):
                seq.id = None
                seq.primers = []
                self.registry.register(sample, seq)

    def _all_dna_samples(self):
        return self.session.Sample.where(
            {"sample_type_id": [d.id for d in self.dna_types]}
        )

    def _all_unregistered_samples(self):
        samples = self._all_dna_samples()
        unregistered = []
        registered = []
        for s in samples:
            e = self.registry.find_in_cache(s)
            if e:
                registered.append(s)
            else:
                unregistered.append(s)
        return unregistered

        # self.info("Resolving sequence for '{}'".format(sample.name))
        # if sample.sample_type_id == self.plasmid_type.id:
        #     seq = self.resolve_plasmid(sample)
        # elif sample.sample_type_id == self.fragment_type.id:
        #     seq = self.resolve_fragment(sample)
        # self.resolved_sequences[sample.id] = seq
        # seq.name = sample.name
        # seq.description = sample.description
        # if 'seq_' not in seq.id:
        #     seq.id = None
        # if hasattr(seq, 'primers') and seq.primers:
        #     seq.primers = []
        # return seq
