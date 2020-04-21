"""
klavinslab

Methods and classes specific to the KlavinsLab
"""

from benchlingapi import Session as BenchlingSession
from pydent import AqSession
import json
from os.path import abspath, dirname, join
from aqbt.aquarium import pydent_utils
from .registry import RegistryConnector, KlavinsLabRegistry
from .registry_crawler import RegistryCrawler
from .linter import Linter
import re
from pydent import models
from typing import List
from aqbt.aquarium.pydent_utils import Constants as C
from pydent import Planner
from primer3plus import Design
from aqbt.bioadapter import convert

here = abspath(dirname(__file__))

with open(join(here, "../secrets/credentials.json"), "r") as f:
    secrets = json.load(f)


class Sessions:
    def __init__(self):
        self._aqproduction = None
        self._aqnursery = None
        self._aqlocal = None
        self._klavinslab_registry_connector = None
        self._klregistry = None
        self._benchling = None

    @property
    def benchling(self):
        if not self._benchling:
            self._benchling = BenchlingSession(secrets["benchling"]["api_key"])
        return self._benchling

    @property
    def aqproduction(self):
        if not self._aqproduction:
            self._aqproduction = AqSession(
                login=secrets["aquarium"]["production"]["login"],
                password=secrets["aquarium"]["production"]["password"],
                aquarium_url=secrets["aquarium"]["production"]["url"],
                name="",
            )  #: production Aquarium session
        return self._aqproduction

    @property
    def aqnursery(self):
        if not self._aqnursery:
            self._aqnursery = AqSession(
                login=secrets["aquarium"]["nursery"]["login"],
                password=secrets["aquarium"]["nursery"]["password"],
                aquarium_url=secrets["aquarium"]["nursery"]["url"],
                name="",
            )  #: production Aquarium session
        return self._aqnursery

    @property
    def aqlocal(self):
        if not self._aqlocal:
            self._aqlocal = AqSession(
                login=secrets["aquarium"]["local"]["login"],
                password=secrets["aquarium"]["local"]["password"],
                aquarium_url=secrets["aquarium"]["local"]["url"],
                name="local",
            )
        return self._aqlocal

    @property
    def klavinslab_registry_connector(self):
        if not self._klavinslab_registry_connector:
            self._klavinslab_registry_connector = RegistryConnector(
                api=self.benchling,
                initials=secrets["benchling_registry"]["initials"],
                schema=secrets["benchling_registry"]["schema"],
                prefix=secrets["benchling_registry"]["prefix"],
                folder_id=secrets["benchling_registry"]["folder_id"],
                registry_id=secrets["benchling_registry"]["id"],
            )  #: klavins lab Benchling registry connector
        return self._klavinslab_registry_connector

    @property
    def klregistry(self):
        if not self._klregistry:
            self._klregistry = KlavinsLabRegistry(
                connector=self.klavinslab_registry_connector,
                aqsession=self.aqproduction,
            )  #: klavins lab Benchling registry
        return self._klregistry


def is_sequencing_file(filename):
    """Check if filename is a sequencing file"""
    if filename.endswith(".ab1"):
        return True
    return False


def parse_seq_filename(filename):
    """Parse aquarium-style sequencing file names."""
    pattern = r"(\d+)-(\w+?)-(\d+)"
    m = re.match(pattern, filename)
    item_id = m.group(1)
    user = m.group(2)
    primer_id = m.group(3)
    return {"item_id": int(item_id), "user": user, "primer_id": int(primer_id)}


def submit_alignments(
    registry: KlavinsLabRegistry, session: AqSession, uploads: List[models.Upload]
):
    """
    Submit alignments to benchling.

    param registry: klregistry
    param session: Aquarium session
    param uploads: list of Upload instances
    """
    groups = {}

    for u in uploads:
        filename = u.download()
        parsed = parse_seq_filename(u.upload_file_name)
        item_id = parsed["item_id"]
        primer_id = parsed["primer_id"]
        user = parsed["user"]
        groups.setdefault(item_id, {})
        if filename in groups[item_id]:
            filename += "_"
        groups[item_id][filename] = {
            "user": user,
            "primer": session.Sample.find(primer_id),
        }

    tasks = {}

    api = registry.connector.api

    for k, v in groups.items():
        print(k)
        item = session.Item.find(k)
        sample = item.sample
        sequence = registry.get_sequence(sample)
        filepaths = [_v for _v in v]
        template = sequence.id

        print(sequence.id)

        task = api.DNAAlignment.submit_alignment(
            algorithm="mafft",
            name="UWBF_{}_{}".format(item.id, sample.name),
            template=sequence.id,  # sequence id of the template
            filepaths=filepaths,
            rawfiles=None,  # only use if you have base64 data handy
        )

        tasks.setdefault(k, list())
        tasks[k].append(task)
    return tasks


def trace_gibson_assembly(sample: models.Sample):
    pass

def make_aq_fragment(
    aqsession,
    template_record,
    name: str,
    project: str,
    description: str,
    template_sample: models.Sample,
    left_overhang: str = "",
    right_overhang: str = "",
):

    # run primer3
    design = Design()
    design.settings.template(str(template_record.seq))
    design.settings.as_cloning_task()
    design.settings.pick_anyway()
    if left_overhang:
        design.settings.left_overhang(str(left_overhang.seq))
    if right_overhang:
        design.settings.right_overhang(str(right_overhang.seq))
    if right_overhang or left_overhang:
        design.settings.use_overhangs()
    primers = design.run()[0][0]

    # check
    assert primers["LEFT"]["location"][0] == 0
    assert primers["RIGHT"]["location"][0] == len(template_record.seq) - 1

    # create record from PCR
    fragment_type = aqsession.SampleType.find_by_name(C.FRAGMENT)
    primer_type = aqsession.SampleType.find_by_name(C.PRIMER)
    record = left_overhang + convert(template_record, to="SeqRecord") + right_overhang

    # make aquarium primers and fragment
    frag = aqsession.Sample.new(
        name=name,
        project=project,
        description=description,
        sample_type=fragment_type,
        properties={
            "Forward Primer": aqsession.Sample.new(
                name=name + "_fwd",
                project=project,
                description="",
                sample_type=primer_type,
                properties={
                    "Anneal Sequence": primers["LEFT"]["SEQUENCE"],
                    "Overhang Sequence": primers["LEFT"]["OVERHANG"],
                    "T Anneal": primers["LEFT"]["TM"] - 3,
                },
            ),
            "Reverse Primer": aqsession.Sample.new(
                name=name + "_rev",
                project=project,
                description="",
                sample_type=primer_type,
                properties={
                    "Anneal Sequence": primers["RIGHT"]["SEQUENCE"],
                    "Overhang Sequence": primers["RIGHT"]["OVERHANG"],
                    "T Anneal": primers["RIGHT"]["TM"] - 3,
                },
            ),
            "Template": template_sample,
            "Length": len(record.seq),
        },
    )
    frag.record = record
    return frag


def submit_alignments_for_plan(session: AqSession, plan_id: int) -> List:
    planner = Planner(session.Plan.find(plan_id))

    uploads = []

    for da in planner.plan.data_associations:
        if da.upload:
            if is_sequencing_file(da.upload.upload_file_name):
                uploads.append(da.upload)

    tasks = submit_alignments(registry, session, uploads)
    return tasks


sessions = Sessions()
