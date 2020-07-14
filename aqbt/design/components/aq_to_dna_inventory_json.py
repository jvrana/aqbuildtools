"""
Inputs:

(1) DNA Registry (e.g. Benchling) + config
(2) Aquarium + config

Outputs:

(1) Inventory JSON for available DNA and their sequences
"""

from aqbt.aquarium.registry import BenchlingRegistryConnector
from pydent import AqSession
from typing import List, Dict, Any, Union
from aqbt.utils import merge_left
from aqbt.bioadapter.conversion import benchling_to_json
import re


def dna_json_from_benchling_connector(
    connector: BenchlingRegistryConnector, limit: int = None
):
    return [benchling_to_json(d) for d in connector.all(limit)]


def available_sample_json_from_aquarium_session(
    session: AqSession,
    sample_type_names: List[str],
    object_type_names: List[str],
    page_size=1000,
    sample_limit: int = -1,
    item_limit: int = -1,
    sample_ids: List[int] = None,
) -> Dict[str, Any]:

    sample_query = {
        "__model__": "Sample",
        "__description__": "get samples",
        "__json__": True,
        "__query__": {
            "sample_type": {"__query__": {"name": sample_type_names}},
            "__options__": {"pageSize": page_size, "limit": sample_limit},
            "__return__": {"field_values": [], "sample_type": [],},
        },
    }
    if sample_ids:
        sample_query["__query__"]["id"] = sample_ids

    samples = session.query(sample_query)

    available_items = session.query(
        {
            "__model__": "Item",
            "__description__": "get available fragment stocks",
            "__json__": True,
            "__query__": {
                "object_type": {"__query__": {"name": object_type_names}},
                "location": {"__not__": "deleted"},
                "sample_id": [entry["id"] for entry in samples],
                "__options__": {"pageSize": page_size, "limit": item_limit},
            },
        }
    )

    samples_with_inventory = merge_left(
        samples,
        available_items,
        "id",
        "sample_id",
        "available_items",
        inplace=True,
        merged_only=True,
    )
    return samples_with_inventory


def sample_ids_from_dna_json(dna_json: List[Dict[str, Any]]):
    pattern = "aqUWBF(?P<aquarium_id>\d+)"
    expected_re_key = "aquarium_id"
    registry_key = "entityRegistryId"
    sample_ids = []
    for entry in dna_json:
        if registry_key in entry:
            registry_id = entry[registry_key]
            matched = re.match(pattern, registry_id)
            if matched:
                sample_ids.append(int(matched.groupdict()[expected_re_key]))
    return sample_ids


# TODO: pattern for converting entityRegistryId to aquarium sample id
def merge_dna_json_and_sample_json(
    dna_json: List[Dict[str, Any]], sample_inventory: List[Dict[str, Any]]
) -> List[Dict[str, Any]]:
    sample_ids = sample_ids_from_dna_json(dna_json)
    dna_json_copy = [dict(d) for d in dna_json]
    for sample_id, entry in zip(sample_ids, dna_json_copy):
        entry["aquarium_sample_id"] = sample_id
    merged = merge_left(
        sample_inventory,
        dna_json_copy,
        "id",
        "aquarium_sample_id",
        "sequence",
        inplace=False,
        merged_only=True,
    )
    return merged


def generate_inventory_json(
    connector: BenchlingRegistryConnector,
    session: AqSession,
    aquarium_dna_types: dict,
    page_size=1000,
    dna_limit: Union[int, None] = None,
    sample_limit: int = -1,
    item_limit: int = -1,
):
    """
    Generate inventory JSON.

    .. code-block:: python

        aquarium_dna_types = {
            'Fragment': {
                'object_types': ['Fragment Stock'],
                'topology': 'linear',
            },
            'Plasmid': {
                'object_types': ['Plasmid Glycerol Stock'],
                'topology': 'cyclic'
            }
        }

    :param connector:
    :param session:
    :param aquarium_dna_types:
    :param page_size:
    :param dna_limit:
    :param sample_limit:
    :param item_limit:
    :return:
    """
    dna_json = dna_json_from_benchling_connector(connector, limit=dna_limit)
    sample_ids = sample_ids_from_dna_json(dna_json)

    # key
    CIRCULAR = "circular"
    LINEAR = "linear"
    TOPOLOGY = "topology"
    OBJECT_TYPES = "object_types"
    BENCHLING_ISCIRCULAR = "isCircular"

    for sample_type_name, types in aquarium_dna_types.items():
        assert OBJECT_TYPES in types
        assert TOPOLOGY in types
        assert types[TOPOLOGY] in [CIRCULAR, LINEAR]

    merged = []
    for sample_type_name, types in aquarium_dna_types.items():
        object_type_names = types[OBJECT_TYPES]
        _sample_json = available_sample_json_from_aquarium_session(
            session,
            sample_type_names=[sample_type_name],
            object_type_names=object_type_names,
            page_size=page_size,
            sample_limit=sample_limit,
            sample_ids=sample_ids,
            item_limit=item_limit,
        )
        _merged = merge_dna_json_and_sample_json(dna_json, _sample_json)

        # correct topologies
        topology = types[TOPOLOGY]

        def _add_sequence_warning(entry, msg):
            warning_key = "__warning__"
            entry.setdefault(warning_key, list())
            entry[warning_key].append(msg)

        for entry in _merged:

            seq = entry["sequence"][0]
            seq_is_circular = seq[BENCHLING_ISCIRCULAR]
            if topology == LINEAR:
                if seq_is_circular:
                    _add_sequence_warning(
                        entry,
                        "DNA sequence '{}' is annotated with `{}=={}` but Aquarium "
                        "sample is indicated as a linear sequence".format(
                            seq["id"], BENCHLING_ISCIRCULAR, seq[BENCHLING_ISCIRCULAR]
                        ),
                    )
                entry["sequence"][0][BENCHLING_ISCIRCULAR] = False
            elif topology == CIRCULAR:
                if not seq_is_circular:
                    _add_sequence_warning(
                        entry,
                        "DNA sequence '{}' is annotated with `{}=={}` but Aquarium "
                        "sample is indicated as a circular sequence".format(
                            seq["id"], BENCHLING_ISCIRCULAR, seq[BENCHLING_ISCIRCULAR]
                        ),
                    )
                entry["sequence"][0][BENCHLING_ISCIRCULAR] = True
        merged += _merged
    return merged
