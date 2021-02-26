import json

import pytest

from aqbt.bioadapter.conversion import json_to_seqrecord
from aqbt.design.components.aq_to_dna_inventory_json import (
    available_sample_json_from_aquarium_session,
)
from aqbt.design.components.aq_to_dna_inventory_json import (
    dna_json_from_benchling_connector,
)
from aqbt.design.components.aq_to_dna_inventory_json import generate_inventory_json
from aqbt.design.components.aq_to_dna_inventory_json import sample_ids_from_dna_json


@pytest.mark.parametrize(
    "x", [(["Fragment"], ["Fragment Stock"]), (["Plasmid"], ["Plasmid Glycerol Stock"])]
)
def test_available_inventory(aquarium, x):

    samples = available_sample_json_from_aquarium_session(
        aquarium, x[0], x[1], sample_limit=100, item_limit=100, page_size=1000
    )
    for s in samples:
        assert "available_items" in s
        assert isinstance(s["available_items"], list)
        assert len(s["available_items"]) > 0


def test_get_all_registered(registry):

    connector = registry.connector
    results = connector.all(None)
    for _ in range(100):
        dna = next(results)
        print(dna)


def test_get_registered_sequences(registry):
    connector = registry.connector
    dna_json = dna_json_from_benchling_connector(connector, limit=10)
    assert dna_json
    assert isinstance(dna_json, list)
    print(dna_json)


# TODO: generate faked DNA JSON
def test_registered_aq_inventory(registry):
    connector = registry.connector
    dna_json = dna_json_from_benchling_connector(connector)
    sample_ids = sample_ids_from_dna_json(dna_json)
    assert sample_ids
    print(sample_ids)


def test_available_inventory(registry, aquarium):
    connector = registry.connector

    result = generate_inventory_json(
        connector,
        aquarium,
        {
            "Fragment": {
                "object_types": ["Fragment Stock"],
                "topology": "linear",
            },
            "Plasmid": {
                "object_types": ["Plasmid Glycerol Stock"],
                "topology": "circular",
            },
        },
        sample_limit=-1,
        item_limit=-1,
        page_size=1000,
        dna_limit=100,
    )
    assert result
    print(json.dumps(result, indent=2))


def test_inventory_to_seq_records(registry, aquarium):
    connector = registry.connector

    result = generate_inventory_json(
        connector,
        aquarium,
        {
            "Fragment": {
                "object_types": ["Fragment Stock"],
                "topology": "linear",
            },
            "Plasmid": {
                "object_types": ["Plasmid Glycerol Stock"],
                "topology": "circular",
            },
        },
        sample_limit=-1,
        item_limit=-1,
        page_size=1000,
        dna_limit=100,
    )

    cyclic_records = [
        json_to_seqrecord(entry["sequence"][0])
        for entry in result
        if entry["sequence"][0]["isCircular"]
    ]
    linear_records = [
        json_to_seqrecord(entry["sequence"][0])
        for entry in result
        if entry["sequence"][0]["isCircular"]
    ]
    assert cyclic_records
    assert linear_records
