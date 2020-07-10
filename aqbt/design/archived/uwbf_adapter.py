import os
import random
from typing import Dict
from typing import List
from typing import Tuple
from typing import Union
from uuid import uuid4

import dasi
import networkx as nx
import pandas as pd
from Bio.SeqRecord import SeqRecord
from dasi.cost import cached_span_cost
from dasi.cost import SpanCost
from lobio import biopython
from lobio.bioadapter import convert
from lobio.klavinslab.pydent_utils import Constants as C
from lobio.klavinslab.registry import KlavinsLabRegistry
from lobio.utils import remove_none_values
from pydent import models
from pydent.models import Sample


DASI_DEFAULT_AQ_TIMEOUT = 120
CYCLIC = "cyclic"
LINEAR = "linear"
IS_DESIGN = "is_final_design"
FLOAT_DECIMAL_PLACES = 2


def create_sample_graph(design: Union[dasi.Design, dasi.LibraryDesign]) -> nx.DiGraph:
    """


    :param design: dasi.Design
    :return:
    """

    def mol_key(mol):
        return (mol.type.name, str(mol.sequence.seq).upper())

    graph = nx.DiGraph()

    for qk, result in design.results.items():
        if result.assemblies:
            assembly = result.assemblies[0]
            reactions = assembly.reactions
            graph.add_node(
                qk, type="design", sequence=design.seqdb[qk], reactions=reactions
            )
            for reaction in reactions:
                molecules = reaction.inputs + reaction.outputs
                for m in molecules:
                    k = mol_key(m)
                    if k in graph:
                        graph.nodes[k]["reactions"].append((qk, reaction, m))
                    else:
                        graph.add_node(
                            k, reactions=[(qk, reaction, m)], type=k[0], sequence=k[1]
                        )
                for mol_out in reaction.outputs:
                    k_out = mol_key(mol_out)
                    graph.add_edge(k_out, qk)
                    for mol_in in reaction.inputs:
                        k_in = mol_key(mol_in)
                        graph.add_edge(k_in, k_out)
        else:
            design.logger.error("No assemblies were found for {}".format(qk))
    return graph


class DesignFactory:

    # TODO: add debug_n_fragments
    def __init__(
        self,
        registry: KlavinsLabRegistry,
        directory: str = ".",
        debug_mode: bool = False,
        debug_n_seqs: int = 100,
        debug_primer_lim: int = 100,
        inventory_df: pd.DataFrame = None,
        ignore_failed_sequences: bool = False,
        span_cost: SpanCost = None,
    ):
        """

        :param registry:
        :param directory:
        :param debug_mode:
        :param debug_n_seqs:
        :param debug_primer_lim:
        :param inventory_df: optional inventory df. If None, will produce a new inventory_df.
        """

        # settings
        self.dir = directory  #: the directory to save
        self.do_load = True  #: if True, will load SpanCost bytes, if available
        self.do_save = True  #: if True, will save SpanCost bytes, if available
        self.debug_n_seqs = (
            debug_n_seqs  #: number of sequences to run while in debug mode
        )
        self.debug_mode = debug_mode  #: run the design in debug mode, generating random Samples, Items, and sequences

        self.registry = registry

        # make inventory df
        if inventory_df is not None:
            self._inventory_df = inventory_df
        elif self.debug_mode:
            self._inventory_df = build_df(
                self.registry,
                dna_limit=debug_n_seqs,
                primer_limit=debug_primer_lim,
                ignore=ignore_failed_sequences,
            )
        else:
            self._inventory_df = build_df(
                self.registry,
                dna_limit=None,
                primer_limit=None,
                ignore=ignore_failed_sequences,
            )

        # default sample types
        self.fragment_type = self.registry.session.SampleType.find_by_name(C.FRAGMENT)
        self.plasmid_type = self.registry.session.SampleType.find_by_name(C.PLASMID)
        self.primer_type = self.registry.session.SampleType.find_by_name(C.PRIMER)

        # load span cost
        print("Loading DASi assembly costs. Please wait.")
        if span_cost is None:
            self.span_cost = cached_span_cost()
        else:
            self.span_cost = span_cost
        print("Assembly costs loaded.")

    @classmethod
    def fake_init(
        cls,
        registry: KlavinsLabRegistry,
        n_plasmids: int = 50,
        n_fragments: int = 50,
        n_primers: int = 50,
        span_cost: SpanCost = None,
    ):
        inv_df = build_fake_inventory_df(registry, n_plasmids, n_fragments, n_primers)
        factory = cls(registry, inventory_df=inv_df, span_cost=span_cost)
        return factory

    @property
    def inventory_df(self) -> pd.DataFrame:
        """
        A pandas.DataFrame that contains inventory information
        for all the samples that are available. the following columns:

        ::

            ALL_COLUMNS = [
                "sample_id",            # the sample id of the Aquarium sample
                "benchling_sequence",   # the Benchling sequence (benchling.DNASequence)
                "record",               # the biopython sequence (Bio.SeqRecord)
                "entity_registry_id",   # the entity_registry_id
                "is_circular",          # whether the topology of the sequence
                "is_available",         # whether there is available inventory for the sample
                "sample",               # the Aquarium Sample instance
                "sample_type",          # the name of the SampleType
                "sequence_hash",        # a SHA1 hash of the sequence for sequence comparison
                "record_uuid",          # the unique ID for the record / sample
            ]

        :return: the dataframe
        """
        return self._inventory_df[self._inventory_df["is_available"] == True]

    def _record_to_row(self, record: SeqRecord) -> pd.Series:
        """Using a Record, find the row in the inventory df that corresponds to
        the LIMS_ID in the record.annotations."""
        sample_id = record.annotations[LIMS_ID]
        return self.inventory_df[self.inventory_df["sample_id"] == sample_id].iloc[0]

    def _record_to_sample(self, record: SeqRecord) -> Sample:
        """Using a Sample, find the row in the inventory df that corresponds to
        the LIMS_ID in the record.annotations.

        :param record:
        :return:
        """
        return self._record_to_row(record)["sample"]

    def _record_seq_to_df(self, record: SeqRecord) -> pd.DataFrame:
        """Using the sha1 hash of the record sequence, find the rows in the
        inventory df that matches the provided record.

        :param record:
        :return:
        """
        hash = seq_sha1(str(record.seq))
        return self.inventory_df[self.inventory_df["sequence_hash"] == hash]

    def get_span_cost(self):
        """Saves the span cost as bytes; reloads when called."""
        path = os.path.join(self.dir, "span_cost.b")
        if self.do_load and os.path.isfile(path):
            with dasi.logger.timeit("INFO", "loading bytes"):
                print("Loading file: {}".format(path))
                span_cost = SpanCost.load(path)
        else:
            span_cost = SpanCost.open()
            if self.do_save:
                with dasi.logger.timeit("INFO", "saving bytes"):
                    print("Saving file: {}".format(path))
                    span_cost.dump(path)
        return span_cost

    def _records(self, sample_type: str):
        df = self.inventory_df
        return list(df[df["sample_type"] == sample_type]["record"])

    @property
    def linear_fragment_records(self):
        return self._records("Fragment")

    @property
    def primer_records(self):
        return self._records("Primer")

    @property
    def cyclic_plasmid_records(self):
        return self._records("Plasmid")

    def _add_materials(
        self, design: Union[dasi.Design, dasi.LibraryDesign], goals: List[SeqRecord]
    ):
        """Add material to the DASi Design instance.

        :param design:
        :param goals:
        :return:
        """
        design.add_fragments(self.linear_fragment_records)
        design.add_primers(self.primer_records)
        design.add_templates(self.cyclic_plasmid_records)
        design.add_queries(goals)

    def _new_library_design(self, goals: List[SeqRecord]) -> dasi.LibraryDesign:
        design = dasi.LibraryDesign(span_cost=self.span_cost)
        self._add_materials(design, goals)
        return design

    def _new_design(self, goals: List[SeqRecord]) -> dasi.Design:
        design = dasi.Design(span_cost=self.span_cost)
        self._add_materials(design, goals)
        return design

    def _goals_from_benchling(self, seq_ids: List[str]) -> List[SeqRecord]:
        records = []
        for seq_id in seq_ids:
            seq = self.registry.connector.api.DNASequence.find(seq_id)
            record = convert(seq, to="SeqRecord")
            if seq.is_circular:
                biopython.make_cyclic([record])
            else:
                biopython.make_linear([record])
            records.append(record)
        return records

    @staticmethod
    def _new_sample_name(qk, mtype):
        """Create a new sample name using the query_key and type_name. Adds a
        unique string value to the end of each name.

        :param qk:
        :param mtype:
        :return:
        """
        return "{}_{}_{}".format(mtype, qk[-5:], str(uuid4())[-4:])

    @staticmethod
    def _random_record_from_library(
        records: List[SeqRecord],
        circular: bool,
        plasmid_size_interval: Tuple[int, int] = (5000, 10000),
        chunk_size_interval: Tuple[int, int] = (100, 3000),
        random_chunk_prob_int: Tuple[float, float] = (0, 0.5),
        random_chunk_size_int: Tuple[int, int] = (100, 1000),
    ):
        """

        :param records:
        :param plasmid_size_interval:
        :param chunk_size_interval:
        :param random_chunk_prob_interval: picks a random probability at whic
        :param circular:
        :return:
        """
        length = random.randint(*plasmid_size_interval)
        new_record = biopython.new_sequence("")

        x = random.uniform(*random_chunk_prob_int)

        while len(new_record.seq) < length:
            if random.random() < x:
                rand_chunk_len = random.randint(*random_chunk_size_int)
                rand_rec = biopython.random_record(rand_chunk_len)
            else:
                rand_chunk_len = random.randint(*chunk_size_interval)
                rand_rec = random.choice(records)
                if rand_chunk_len > len(rand_rec.seq):
                    continue
                    # rand_chunk_len = len(rand_rec.seq)
            i = random.randint(0, len(rand_rec.seq) - rand_chunk_len)
            rand_chunk = rand_rec[i : i + rand_chunk_len]
            if random.random() < 0.5:
                rand_chunk = rand_chunk.reverse_complement()
            new_record += rand_chunk

        new_record = new_record[:length]
        if circular:
            biopython.make_cyclic([new_record])
        else:
            biopython.make_linear([new_record])
        return new_record

    def design_fake(
        self,
        n_fake_designs: int = 3,
        plasmid_size_interval: Tuple[int, int] = None,
        chunk_size_interval: Tuple[int, int] = None,
        random_chunk_prob_int: Tuple[float, float] = None,
    ):
        """Produce faked records from the inventory_df."""
        kwargs = dict(
            plasmid_size_interval=plasmid_size_interval,
            chunk_size_interval=chunk_size_interval,
            random_chunk_prob_int=random_chunk_prob_int,
        )
        remove_none_values(kwargs)
        records = []
        library_records = list(self.inventory_df["record"])
        for _ in range(n_fake_designs):
            records.append(
                self._random_record_from_library(
                    library_records, circular=True, **kwargs
                )
            )
        return self.design_from_records(records)

    def design_from_records(self, records: List[SeqRecord]):
        design = self._new_library_design(goals=records)
        design.compile_library(n_jobs=1)
        design.optimize_library()

        graph = create_sample_graph(design)
        all_samples, design_dict = self.design_to_samples_json(graph)
        return design, all_samples, design_dict, graph

    def design_from_benchling(
        self, seq_ids: List[str]
    ) -> Tuple[
        Union[dasi.Design, dasi.LibraryDesign], List[models.Sample], Dict, nx.DiGraph
    ]:

        recs = self._goals_from_benchling(seq_ids)
        return self.design_from_records(recs)

    def design_to_samples_json(
        self, graph: nx.DiGraph
    ) -> Tuple[List[models.Sample], Dict]:
        """

        :param graph:
        :return:
        """
        for node, ndata in graph.nodes(data=True):
            ndata["sample"] = None

        for node, ndata in graph.nodes(data=True):
            mtype = ndata["type"]
            ndata["sample"] = None
            if mtype in ["TEMPLATE"]:
                self.resolve_template(ndata)
            elif mtype in ["PRIMER"]:
                self.resolve_primer(ndata)
            elif mtype in ["PRE-MADE DNA FRAGMENT"]:
                self.resolve_premade_fragment(ndata)
            elif mtype in ["GAP", "SHARED_SYNTHESIZED_FRAGMENT"]:
                self.resolve_gap(ndata)

        for node, ndata in graph.nodes(data=True):
            mtype = ndata["type"]
            if "PCR" in mtype:
                self.resolve_pcr(ndata, node, graph)
            elif mtype == "design":
                self.resolve_design(ndata)

        all_samples = []
        sample_node_to_index = {}
        for i, (node, ndata) in enumerate(graph.nodes(data=True)):
            all_samples.append(ndata["sample"])
            sample_node_to_index[node] = i

        design_dict = {}

        for node, ndata in graph.nodes(data=True):
            i1 = sample_node_to_index[node]
            design_dict[i1] = {"sample": ndata["sample"], "inputs": list()}
            if IS_DESIGN in ndata:
                for n in graph.predecessors(node):
                    pndata = graph.nodes[n]
                    design_dict[i1]["inputs"].append(pndata["sample"])

        return all_samples, design_dict

    def resolve_template(self, ndata):
        sorted_reactions = sorted(
            ndata["reactions"],
            key=lambda x: x[2].alignment_group.query_region.cyclic,
            reverse=True,
        )
        selected = sorted_reactions[0]
        molecule = selected[2]
        record = molecule.sequence
        aqsample = self._record_to_sample(record)
        ndata["sample"] = aqsample

    def new_primer(
        self, anneal: str, overhang: str, t_anneal: float, name: str
    ) -> Sample:
        """Makes a new Aquarium primer.

        :param anneal: anneal sequence
        :param overhang: overhang sequence
        :param t_anneal: annealing temperature
        :param name: name of primer
        :return: Sample
        """
        return self.primer_type.new_sample(
            name=name,
            description="",
            project="SD2",
            properties={
                "Anneal Sequence": anneal,
                "Overhang Sequence": overhang,
                "T Anneal": round(t_anneal, FLOAT_DECIMAL_PLACES),
            },
        )

    def resolve_primer(self, ndata: Dict):
        """Finds (or creates a new) Aquarium primer.

        :param ndata: node data
        :return:
        """
        selected = ndata["reactions"][0]
        molecule = selected[2]
        record = molecule.sequence
        try:
            aqprimer = self._record_seq_to_df(record).iloc[0]["sample"]
            ndata["sample"] = aqprimer
        except IndexError:
            meta = molecule.metadata
            qk = ndata["reactions"][0][0]
            primer_name = self._new_sample_name(qk, "primer")
            ndata["sample"] = self.new_primer(
                meta["SEQUENCE"],
                meta["OVERHANG"],
                float("%.2f" % round(meta["TM"] - 3, 2)),
                primer_name,
            )

    def resolve_premade_fragment(self, ndata: Dict):
        """Finds existing Aquarium fragment.

        :param ndata: node data
        :return: None
        """
        selected = ndata["reactions"][0]
        molecule = selected[2]
        record = molecule.sequence
        aqsample = self._record_to_sample(record)
        ndata["sample"] = aqsample

    def resolve_pcr(self, ndata: Dict, node: str, graph: nx.DiGraph):
        """Creates a new Aquarium sample for a PCR reaction.

        :param ndata: node data
        :param node: node identifier
        :param graph: directed graph
        :return: None
        """
        selected = ndata["reactions"][0]
        molecule = selected[2]
        qk = selected[0]
        record = molecule.sequence

        preds = graph.predecessors(node)
        templates = []
        primers = []
        for n in preds:
            pred_ndata = graph.nodes[n]
            if pred_ndata["type"] == "TEMPLATE":
                templates.append(pred_ndata["sample"])
            elif pred_ndata["type"] == "PRIMER":
                primers.append(pred_ndata["sample"])
        assert len(primers) == 2
        assert len(templates) >= 1

        ndata["sample"] = self.fragment_type.new_sample(
            name=self._new_sample_name(qk, "fragment"),
            description="",
            project="",
            properties={
                "Sequence": str(record.seq),
                "Length": len(record.seq),
                "Template": templates[0],
                "Forward Primer": primers[0],
                "Reverse Primer": primers[1],
            },
        )

    def resolve_gap(self, ndata: Dict):
        """Creates a new Aquarium sample for a 'gap' segment.

        :param ndata: node data
        :return: None
        """
        selected = ndata["reactions"][0]
        qk = selected[0]
        molecule = selected[2]
        record = molecule.sequence
        ndata["sample"] = self.fragment_type.new_sample(
            name=self._new_sample_name(qk, "gblock"),
            project="",
            description="",
            properties={"Sequence": str(record.seq), "Length": len(record.seq)},
        )

    def resolve_design(self, ndata: Dict):
        """Creates a new Aquarium sample for the goal plasmid.

        :param ndata: node data
        :return: None
        """
        ndata[IS_DESIGN] = True
        ndata["sample"] = self.plasmid_type.new_sample(
            name=ndata["sequence"].name,
            description="",
            project="",
            properties={
                "Sequence": ndata["sequence"].id,
                "Length": len(ndata["sequence"].seq),
            },
        )
