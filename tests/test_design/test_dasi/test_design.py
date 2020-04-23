from os.path import abspath
from os.path import dirname

import pytest
from dasi import Design
from dasi.utils.biopython import random_record_from_library

from aqbt.design.dasi.dfs import KlavinsLabDnaDb

# TODO: replace DesignFactory
# from aqbt.design.dasi.uwbf_adapter import DesignFactory


here = abspath(dirname(__file__))


def test_build_fake_df(registry):
    registry = registry.copy()
    db = KlavinsLabDnaDb(registry)
    db.build_fake_inventory_df(50, 50, 50)


def test_small_design(registry):
    registry = registry.copy()
    db = KlavinsLabDnaDb(registry)
    db.build(dna_limit=50, primer_limit=50)

    df = db.df

    def _get_records(st):
        return list(df[df["sample_type"] == st]["record"])

    fragments = _get_records("Fragment")
    plasmids = _get_records("Plasmid")
    primers = _get_records("Primer")

    designs = [random_record_from_library(fragments + plasmids, circular=True)]
    design = Design()
    design.add_materials(
        primers=primers, fragments=fragments, templates=plasmids, queries=designs
    )
    design.run(n_paths=1)
    print(design.out())


def test_large_design(registry):
    registry = registry.copy()
    db = KlavinsLabDnaDb(registry)
    db.build(dna_limit=None, primer_limit=None)

    df = db.df

    def _get_records(st):
        return list(df[df["sample_type"] == st]["record"])

    fragments = _get_records("Fragment")
    plasmids = _get_records("Plasmid")
    primers = _get_records("Primer")

    designs = [
        random_record_from_library(fragments + plasmids, circular=True)
        for _ in range(1)
    ]
    design = Design()
    design.add_materials(
        primers=primers, fragments=fragments, templates=plasmids, queries=designs
    )
    design.run(n_paths=1)
    print(design.out())


def test_dasi(registry):
    registry = registry.copy()
    db = KlavinsLabDnaDb(registry)
    db.build(dna_limit=50, primer_limit=50)


########################
# Test designs
########################

# TODO: replace DesignFactory tests
#
# def test_dasi_design_debug_mode_init(registry, cached_span_cost):
#     registry = registry.copy()
#     factory = DesignFactory(
#         registry,
#         debug_mode=True,
#         debug_n_seqs=50,
#         debug_primer_lim=50,
#         span_cost=cached_span_cost,
#     )
#     assert factory


# def test_dasi_design_production_init(registry, cached_span_cost):
#     registry = registry.copy()
#     registry.session.set_timeout(120)
#     registry.benchling.set_timeout(60)
#     factory = DesignFactory(registry, span_cost=cached_span_cost)
#     assert factory


# # TODO: move the random library code to DASi for testing
# @pytest.mark.parametrize(
#     "random_chunk_prob_int", [(0, 0), (0.1, 0.1), (0.3, 0.3), (0.5, 0.5)]
# )
# @pytest.mark.parametrize("n_builds", [1, 3])
# def test_fake_design(registry, random_chunk_prob_int, n_builds, cached_span_cost):
#     registry = registry.copy()
#     factory = DesignFactory.fake_init(registry, 50, 50, 50, span_cost=cached_span_cost)
#     design, all_samples, design_dict, graph = factory.design_fake(
#         n_builds, random_chunk_prob_int=random_chunk_prob_int
#     )
#     print(design_dict)
#     df1, df2, design_json = design.to_df()
#     for _, d in design_dict.items():
#         print(d["sample"].sample_type.name)
#         print(d["sample"].id)
#     print(df1)
#     print(df2)
#

# def test_design(registry, cached_span_cost):
#     registry = registry.copy()
#     factory = DesignFactory(registry, debug_mode=True, span_cost=cached_span_cost)
#
#     design, all_samples, design_dict, graph = factory.design_from_benchling(
#         ["seq_1i2wmZhy", "seq_FyA4Jvjt"]
#     )
#     for _, d in design_dict.items():
#         print(d["sample"].sample_type.name)
#         print(d["sample"].id)
#     df1, df2, design_json = design.to_df()
#     print(df1)
#     print(df2)
