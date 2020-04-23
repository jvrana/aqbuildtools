from typing import List
from typing import Tuple

import networkx as nx
from pydent import AqSession
from pydent import Planner
from pydent.models import FieldType
from pydent.models import FieldValue
from pydent.models import Operation
from pydent.models import Sample


class Constants:

    # SampleType names
    YEAST = "Yeast Strain"
    PLASMID = "Plasmid"
    FRAGMENT = "Fragment"
    PRIMER = "Primer"

    # ObjectTypes
    PLASMID_GLYCEROL_STOCK = "Plasmid Glycerol Stock"
    FRAGMENT_STOCK = "Fragment Stock"

    # FieldValue values
    DIPLOID = "Diploid"
    MATa = "MATa"
    MATalpha = "MATalpha"

    # FieldType names
    MATING_TYPE = "Mating Type"
    HAPLOIDS = "Haploids"
    INTEGRANT = "Integrant"
    PARENT = "Parent"

    # Other
    DELETED = "deleted"


def io_using_samples(
    session: AqSession, samples: List[Sample], role: str = None
) -> List[FieldValue]:
    query = {"parent_class": "Operation", "child_sample_id": [s.id for s in samples]}
    if role:
        query["role"] = role
    fvs = session.FieldValue.where(query)
    return fvs


def inputs_using_samples(session: AqSession, samples: List[Sample]) -> List[FieldValue]:
    return io_using_samples(session, samples, role="input")


def outputs_using_samples(
    session: AqSession, samples: List[Sample]
) -> List[FieldValue]:
    return io_using_samples(session, samples, role="output")


def ops_using_samples(session, samples, role=None, status=None, group_by_sample=False):
    fvs = io_using_samples(session, samples, role=role)
    query = {"id": [fv.parent_id for fv in fvs]}

    if status:
        query["status"] = status
    ops = session.Operation.where(query)

    if group_by_sample:
        sample_to_ops = {}
        for fv in fvs:
            sample = session.Sample.find(fv.child_sample_id)
            op = session.Operation.find(fv.parent_id)
            sample_to_ops.setdefault(sample, list())
            sample_to_ops[sample].append(op)
        return sample_to_ops
    else:
        return ops


def is_keystone_op(op: Operation) -> bool:
    input_sample_ids = [fv.child_sample_id for fv in op.inputs]
    output_sample_ids = [fv.child_sample_id for fv in op.outputs]
    for sid in output_sample_ids:
        if sid not in input_sample_ids:
            return True
    return False


def pages(interface, query={}, page_size=50, num_pages=2, first=False):
    for page in range(num_pages):
        if first:
            yield interface.first(
                page_size, query=query, opts={"offset": page_size * page}
            )
        else:
            yield interface.last(
                page_size, query=query, opts={"offset": page_size * page}
            )


def sample_pages(
    session: AqSession,
    query: dict = None,
    page_size: int = 50,
    num_pages: int = 2,
    first: bool = False,
):
    return pages(
        session.Sample,
        query=query,
        page_size=page_size,
        num_pages=num_pages,
        first=first,
    )


class Status:
    planning = "planning"
    waiting = "waiting"
    primed = "primed"
    pending = "pending"
    delayed = "delayed"
    done = "done"
    error = "error"
    scheduled = "scheduled"

    statuses = [error, waiting, delayed, pending, scheduled, planning, primed, done]
    status_level = list(range(len(statuses)))[::-1]
    level_to_status = dict(zip(status_level, statuses))
    status_to_level = dict(zip(statuses, status_level))

    @staticmethod
    def sample_status(session: AqSession, samples: List[Sample]) -> dict:
        status_dict = {}
        with session.with_cache(using_models=True, timeout=60) as sess:
            sample_to_ops = ops_using_samples(
                sess, samples, role="output", group_by_sample=True
            )
            for sample, ops in sample_to_ops.items():
                level = [Status.status_to_level[op.status] for op in ops][0]
                status = Status.level_to_status[level]
                status_dict[sample] = {
                    "status": status,
                    "plans": sess.browser.get(ops, "plans"),
                }
        return status_dict

    @staticmethod
    def production_status(session: AqSession, samples: List[Sample]) -> dict:
        """Get the per plan sample production status."""
        with session.with_cache() as sess:
            ops_by_sample = ops_using_samples(
                sess, samples, role="output", group_by_sample=True
            )
            sample_production_status = {}
            for sample, ops in ops_by_sample.items():
                keystone_ops = [op for op in ops if is_keystone_op(op)]
                plans = sess.browser.get(keystone_ops, "plans")
                for plan in plans:
                    sess.browser.get(plan.operations, "operation_type")
                    planner = Planner(plan)
                    graph = planner.graph
                    plan_ops = [
                        op
                        for op in keystone_ops
                        if graph.ops_to_nodes([op])[0] in graph.nxgraph
                    ]

                    nodes = []
                    for op in plan_ops:
                        nxgraph = nx.DiGraph()
                        nxgraph.add_edges_from(graph.nxgraph.edges(data=False))
                        nxgraph.add_nodes_from(graph.nxgraph.nodes())

                        tree = nx.bfs_tree(
                            nxgraph.reverse(), graph.ops_to_nodes([op])[0]
                        )
                        nodes += list(tree.nodes())

                    nodes = [
                        n for n in nx.topological_sort(graph.nxgraph) if n in nodes
                    ]

                    pred_ops = graph.nodes_to_ops(nodes)

                    #             pred_ops = graph.ops_to_predecessors(plan_ops)
                    sample_production_status.setdefault(sample, {})
                    status_level = [
                        Status.status_to_level[op.status] for op in pred_ops + plan_ops
                    ]
                    overall_status = Status.level_to_status[max(status_level)]
                    grouped_by_status = {}
                    for op in pred_ops + plan_ops:
                        grouped_by_status.setdefault(op.status, list())
                        grouped_by_status[op.status].append(op)

                    sample_production_status[sample][plan.id] = {
                        "planner": planner,
                        "dependent_ops": pred_ops,
                        "ops": plan_ops,
                        "status": overall_status,
                        "culprits": grouped_by_status[overall_status],
                    }
        return sample_production_status


def _get_choices(ft: FieldType):
    return ft.choices.split(",")


def _validate_session(session: AqSession):
    _validate_yeast_type(session)


def _get_session_from_sample(sample: Sample):
    session = sample.session
    _validate_session(session)
    return session


def _validate_yeast_type(session: AqSession):
    # validate Yeast type
    yeast = session.SampleType.find_by_name(Constants.YEAST)
    ft = yeast.field_type(Constants.MATING_TYPE)
    assert Constants.DIPLOID in _get_choices(ft)

    assert yeast.field_type(Constants.HAPLOIDS)


def is_yeast(sample: Sample) -> bool:
    session = _get_session_from_sample(sample)
    return sample.sample_type_id == session.SampleType.find_by_name(Constants.YEAST).id


def is_diploid(sample: Sample) -> bool:
    assert is_yeast(sample)
    session = _get_session_from_sample(sample)
    _validate_session(session)
    return sample.properties[Constants.MATING_TYPE] == Constants.DIPLOID


def is_mat_alpha(sample: Sample) -> bool:
    assert is_yeast(sample)
    session = _get_session_from_sample(sample)
    _validate_session(session)
    return sample.properties[Constants.MATING_TYPE] == Constants.MATalpha


def is_mat_a(sample: Sample) -> bool:
    assert is_yeast(sample)
    session = _get_session_from_sample(sample)
    _validate_session(session)
    return sample.properties[Constants.MATING_TYPE] == Constants.MATa


def is_haploid(sample: Sample) -> bool:
    assert is_yeast(sample)
    session = _get_session_from_sample(sample)
    _validate_session(session)
    return sample.properties[Constants.MATING_TYPE] in [
        Constants.MATa,
        Constants.MATalpha,
    ]


def get_haploids(sample: Sample) -> bool:
    assert is_yeast(sample)
    return sample.properties.get(Constants.HAPLOIDS, [])


def get_integrant(sample: Sample) -> Sample:
    assert is_yeast(sample)
    return sample.properties[Constants.INTEGRANT]


def get_parent(sample: Sample) -> Sample:
    assert is_yeast(sample)
    return sample.properties[Constants.PARENT]


def _a_alpha_or_neither(
    samples: List[Sample],
) -> Tuple[List[Sample], List[Sample], List[Sample]]:
    a, alpha, neither = [], [], []
    for s in samples:
        if is_mat_a(s):
            a.append(s)
        elif is_mat_alpha(s):
            alpha.append(s)
        else:
            neither.append(s)
    return a, alpha, neither


def _get_haploid_by_mating_types(
    sample: Sample,
) -> Tuple[List[Sample], List[Sample], List[Sample]]:
    assert is_yeast(sample)
    if not is_diploid(sample):
        return None
    else:
        haploids = get_haploids(sample)
        mata, matalpha, neither = _a_alpha_or_neither(haploids)
        return mata, matalpha, neither


def get_haploid_by_mating_types(sample: Sample) -> List[Sample]:
    mata, matalpha, neither = _get_haploid_by_mating_types(sample)
    return {Constants.MATa: mata, Constants.MATalpha: matalpha, "Neither": neither}
