from copy import deepcopy
from os.path import abspath
from os.path import dirname
from os.path import join
from typing import Any
from typing import Dict
from typing import List
from typing import Tuple

import networkx as nx
from Bio.SeqRecord import SeqRecord
from pydent.models import Sample

from aqbt.aquarium import pydent_utils
from aqbt.aquarium.registry import LabDNARegistry
from aqbt.aquarium.registry import RetrievalPriorities
from aqbt.bioadapter import convert
from aqbt.biopython import make_linear
from aqbt.biopython import to_gff
from aqbt.genome import _generate_integration_graph
from aqbt.genome import _get_integration_fragments
from aqbt.genome import GenomeIntegrator
from aqbt.genome import IntegrationError
from aqbt.logger import logger


# TODO: generic register interface? e.g. SBOL
def integration(
    registry: LabDNARegistry,
    sample: Sample,
    genome_dictionary: Dict[int, List[SeqRecord]],
    restriction_enzyme: str,
    max_integration_overhang: int = 100,
) -> Dict[str, Dict[str, Any]]:
    """

    :param registry:
    :param sample:
    :param genome_dictionary:
    :param restriction_enzyme:
    :param max_integration_overhang:
    :return:
    """
    session = registry.session

    integrator = GenomeIntegrator(max_integration_overhang=max_integration_overhang)

    with session.with_cache(timeout=60) as sess:
        new_graph = _generate_integration_graph(sess, sample)
        yeast_type = sess.SampleType.find_by_name("Yeast Strain")
        nodes = nx.topological_sort(new_graph)

        for _, sample_id in nodes:
            records = genome_dictionary.get(sample_id, {"genome": None})["genome"]
            yeast = sess.Sample.find(sample_id)
            if yeast.sample_type_id == yeast_type.id:
                if not records:
                    plasmid = yeast.properties["Integrant"]
                    parent = yeast.properties["Parent"]
                    parent_records = genome_dictionary.get(parent.id, {"genome": None})[
                        "genome"
                    ]
                    if not parent_records:
                        raise ValueError(
                            "Cannot perform integration. {} missing records".format(
                                parent.id
                            )
                        )
                    plasmid_seq = registry.get_sequence(
                        plasmid,
                        priority=[
                            RetrievalPriorities.CACHE,
                            RetrievalPriorities.REGISTRY,
                            RetrievalPriorities.URL,
                        ],
                    )
                    err_msg = "Could not find an integration event for {}: {}".format(
                        yeast.id, yeast.name
                    )
                    if not plasmid_seq:
                        raise IntegrationError(
                            err_msg + "\nNo sequence found for "
                            "{}: {}".format(plasmid.id, plasmid.name)
                        )
                    print(plasmid_seq.web_url)

                    plasmid_record = convert(plasmid_seq, to="SeqRecord")
                    integration_records = _get_integration_fragments(
                        plasmid_record, restriction_enzyme
                    )

                    # make everything linear
                    make_linear(integration_records)
                    make_linear(parent_records)

                    result = integrator.integrate(parent_records, integration_records)
                    if result["errors"]:
                        raise IntegrationError(
                            "Could not perform integration:\n{}".format(
                                result["errors"]
                            )
                        )
                    genome_dictionary[yeast.id] = {
                        "genome": result["genome"],
                        "integrations": genome_dictionary.get(yeast.id, {}).get(
                            "integrations", [result["integration"]]
                        ),
                    }
    return genome_dictionary


def mating(
    mat_a_records: List[SeqRecord], mat_alpha_records: List[SeqRecord]
) -> List[SeqRecord]:
    """Make a new list of 'mated' records by combining the two list of
    records."""
    records = []

    for r in mat_a_records:
        r = deepcopy(r)
        r.id = "MatA_" + str(r.id)
        r.name = "MatA_" + r.name
        records.append(r)

    for r in mat_alpha_records:
        r = deepcopy(r)
        r.id = "MatAlpha_" + str(r.id)
        r.name = "MatAlpha_" + r.name
        records.append(r)
    return records


def aq_to_gff(
    registry: LabDNARegistry,
    yeast_strain: Sample,
    starting_dictionary: Dict[int, List[SeqRecord]],
    restriction_enzyme="PmeI",
    out_handle: str = None,
    out_dir: str = None,
    do_save: bool = True,
) -> Tuple[str, Dict[int, List[SeqRecord]]]:
    """Produce a GFF from an Aquarium sample.

    :param registry: klavins lab registry
    :param yeast_strain: yeast strain
    :param starting_dictionary: dictionary of aq_sample_ids to list of Records (for
        base strains like W303 or CEN.PK
    :param restriction_enzyme: which restriction enzyme to use to linearize the plasmids (default: PmeI)
    :param out_handle: output filename
    :param out_dir: output directory (optional)
    :return:
    """
    # make filepath
    if out_handle is None:
        out_handle = "UWBF_{}.gff".format(yeast_strain.id)
    if out_dir is None:
        out_dir = dirname(abspath("."))
    filepath = join(out_dir, out_handle)

    trace = []

    # make records
    if pydent_utils.is_diploid(yeast_strain):
        types = pydent_utils.get_haploid_by_mating_types(yeast_strain)
        mata = types[pydent_utils.Constants.MATa]
        matalpha = types[pydent_utils.Constants.MATalpha]
        neither = types["Neither"]
        if neither:
            print("Found haploids that are neither mat a or mat alpha")
            for s in neither:
                print(s.name, s.properties["Mating Type"])
        assert not neither
        assert len(mata) == 1
        assert len(matalpha) == 1
        mata = mata[0]
        matalpha = matalpha[0]
        mat_a_trace = integration(
            registry,
            mata,
            genome_dictionary=starting_dictionary,
            restriction_enzyme=restriction_enzyme,
        )[mata.id]
        mat_a_records = mat_a_trace["genome"]
        mat_alpha_trace = integration(
            registry,
            matalpha,
            genome_dictionary=starting_dictionary,
            restriction_enzyme=restriction_enzyme,
        )[matalpha.id]
        mat_alpha_records = mat_alpha_trace["genome"]
        chr_records = mating(
            mat_a_records=mat_a_records, mat_alpha_records=mat_alpha_records
        )
        trace += [mat_a_trace, mat_alpha_trace]
    else:
        integration_trace = integration(
            registry,
            yeast_strain,
            genome_dictionary=starting_dictionary,
            restriction_enzyme=restriction_enzyme,
        )[yeast_strain.id]
        chr_records = integration_trace["genome"]
        trace.append(integration_trace)
    # save to GFF
    if do_save:
        logger.debug("Saving records to '{}'".format(filepath))
        path = to_gff(records=chr_records, out_path=filepath, in_place=False)
    else:
        path = ""
    return path, starting_dictionary, trace
