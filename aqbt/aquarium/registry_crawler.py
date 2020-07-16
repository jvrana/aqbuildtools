import networkx as nx
from pydent.base import ModelBase
from pydent.models import Sample

from aqbt.aquarium.linter import Linter
from aqbt.aquarium.pydent_utils import Constants as C
from aqbt.aquarium.pydent_utils import is_keystone_op
from aqbt.logger import logger


class RegistryCrawler:
    OP_SAMPLES_KEY = "op_samples"
    OP_KEY = "operation"

    def __init__(self, session, registry):
        self.registry = registry.copy()
        self.registry.using_cache = True
        self.session = session.with_cache(timeout=60, using_models=True)
        self.dna_types = self.session.SampleType.where({"name": [C.PRIMER, C.FRAGMENT]})
        self.logger = logger(self)

    @property
    def valid_dna_sample_type_ids(self):
        return [d.id for d in self.dna_types]

    def dna_like(self, sample):
        return (
            isinstance(sample, Sample)
            and sample.sample_type_id in self.valid_dna_sample_type_ids
        )

    def only_dna_like(self, models):
        return [m for m in models if self.dna_like(m)]

    def get_models(self, model):
        for fv in model.field_values:
            if fv.sample and issubclass(type(fv.sample), ModelBase):
                yield fv.sample, {"field_value_name": fv.name}
        if hasattr(model, self.OP_SAMPLES_KEY):
            op_samples = getattr(model, self.OP_SAMPLES_KEY)
            for op_sample in op_samples:
                yield op_sample, {"from_operation": getattr(model, self.OP_KEY)}

    def op_selector(self, ops):
        ops = [op for op in ops if op.status != "planning"]
        if len(ops) == 1:
            return ops[0]
        elif not ops:
            return None
        else:
            return None

    def cache_plasmids(self, plasmids):
        fvs = self.session.FieldValue.where(
            {
                "parent_class": "Operation",
                "child_sample_id": [p.id for p in plasmids],
                "role": "output",
            }
        )
        ops = self.session.browser.get(
            fvs, {"operation": {"operation_type": {}, "field_values": "sample"}}
        )["operation"]

        grouped = {}
        for op in ops:
            if is_keystone_op(op):
                fvs = [
                    fv
                    for fv in op.outputs
                    if fv.child_sample_id in [p.id for p in plasmids]
                ]
                for fv in fvs:
                    grouped.setdefault(fv.child_sample_id, list())
                    grouped[fv.child_sample_id].append(op)

        selected_grouped = {}
        for sample_id, ops in grouped.items():
            selected_grouped[sample_id] = self.op_selector(ops)

        for sample_id, op in selected_grouped.items():
            if op:
                op_samples = self.only_dna_like([fv.sample for fv in op.inputs])
                sample = self.session.Sample.find(sample_id)
                setattr(sample, self.OP_SAMPLES_KEY, op_samples)
                setattr(sample, self.OP_KEY, op)

    def cache_func(self, models):

        plasmids = [m for m in models if m.sample_type.name == "Plasmid"]
        self.cache_plasmids(plasmids)
        self.session.browser.get(models, {"field_values": "sample", "sample_type": {}})

    def sample_network(self, samples):
        """Build a DAG of :class:`Samples <pydent.models.Sample>` from their.

        :class:`FieldValues <pydent.models.FieldValue>`.

        .. versionadded:: 0.1.5a7
            method added

        :param samples: list of samples
        :param reverse: whether to reverse the edges of the final graph
        :param g: the graph
        :return:
        """

        return self.session.browser.relationship_network(
            samples, get_models=self.get_models, cache_func=self.cache_func
        )

    def crawler_graph_init(self, samples):
        g = self.sample_network(samples)
        for n, ndata in g.nodes(data=True):
            ndata["errors"] = []
            ndata["sequence"] = None
            ndata["is_registered"] = False
        return g

    def crawler_graph_check_registered(self, g):
        for n, ndata in g.nodes(data=True):
            sample = self.session.Sample.find(n[1])
            is_registered = False
            found_seq = self.registry.find_in_cache(sample)
            if found_seq:
                ndata["is_registered"] = True
                ndata["sequence"] = found_seq
            ndata["is_registered"] = is_registered

    def crawler_graph_update_sequences(self, g):
        for n, ndata in g.nodes(data=True):
            sample = self.session.Sample.find(n[1])
            if not ndata["is_registered"]:
                seq = self.registry.get_sequence(sample)
                if seq:
                    ndata["sequence"] = seq

    def crawler_graph_check_errors(self, g):
        for n in nx.topological_sort(g):
            sample = self.session.Sample.find(n[1])
            ndata = g.nodes[n]

            if sample.sample_type.name == "Primer":
                continue
            elif sample.sample_type.name not in [C.PLASMID, C.FRAGMENT]:
                continue

            if not ndata["sequence"]:
                parents = list(g.predecessors(n))

                if not parents:
                    ndata["errors"].append("{} has no predecessors".format(n))
                else:
                    parent_errors = []
                    for parent in parents:
                        pdata = g.nodes[parent]
                        parent_errors += pdata["errors"]
                    if parent_errors:
                        ndata["errors"].append("{} parents have errors".format(n))
                    elif sample.sample_type.name == "Plasmid":
                        pass
                    elif sample.sample_type.name == C.FRAGMENT:
                        if not len(parents) >= 3:
                            ndata["errors"].append("fragment is missing parents")

    def crawler_register_sequences(self, g):
        for n, ndata in g.nodes(data=True):
            sample = self.session.Sample.find(n[1])
            if not ndata["is_registered"] and ndata["sequence"]:
                try:
                    seq = self.registry.register(
                        sample, ndata["sequence"], overwrite=False, do_raise=False
                    )
                    if seq:
                        ndata["sequence"] = seq
                except Exception as e:
                    ndata["errors"].append(str(e))

    def crawler_graph_build_sequences(self, g):
        for n in nx.topological_sort(g):
            ndata = g.nodes[n]
            print(ndata)
            sample = self.session.Sample.find(n[1])
            if not ndata["sequence"] and not ndata["errors"]:
                print("building sequence for {}".format(n))
                parents = g.predecessors(n)
                if sample.sample_type.name == "Plasmid":
                    sequences = [g.nodes[p]["sequence"] for p in parents]
                    print(sequences)
                elif sample.sample_type.name == C.FRAGMENT:
                    linter = Linter()
                    products = linter.lint_fragment(self.registry, sample)
                    linter.report()
                    if linter.errors:
                        ndata["errors"] += linter.errors
                        for e in linter.errors:
                            self.logger.error(e)
                    elif len(products) == 1:
                        product = products[0]
                        ndata["sequence"] = self.registry.connector.convert(
                            product, to="DNASequence"
                        )
                    elif len(products) > 1:
                        self.logger.error(
                            "More than one product found for {}".format(n)
                        )
                        ndata["errors"].append("More than one product")

    # TODO: auto_associate data for where it was derived from
    # TODO: specify benchling registry id
    # TODO: auto register DNA
    def crawl(self, samples):
        self.logger.info("Initializing graph for {} samples...".format(len(samples)))
        g = self.crawler_graph_init(samples)
        self.logger.info(
            "Graph initialized with {} nodes and {} edges".format(
                g.number_of_nodes(), g.number_of_edges()
            )
        )

        self.logger.info("Checking registered DNASequences...")
        self.crawler_graph_check_registered(g)
        self.logger.info(
            "Found {} registered sequences".format(
                len([n for n, ndata in g.nodes(data=True) if ndata["is_registered"]])
            )
        )

        self.logger.info("Updating sequences from sharelinks or sample properties...")
        self.crawler_graph_update_sequences(g)

        self.logger.info("Checking graph for errors...")
        self.crawler_graph_check_errors(g)

        self.logger.info("Building sequences...")
        self.crawler_graph_build_sequences(g)

        self.logger.info("Registering sequences...")
        self.crawler_register_sequences(g)
