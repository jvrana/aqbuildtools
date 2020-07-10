"""faker.py."""
import random
from typing import Tuple
from uuid import uuid4

import primer3
import primer3plus
from pydent import AqSession
from pydent.models import Sample
from pydent.models import SampleType

from aqbt import biopython
from aqbt.aquarium.pydent_utils import Constants

# TODO: config should be used here
class FakeSampleGenerator:
    def __init__(self, session: AqSession):
        self.session = session
        self.plasmid_type = session.SampleType.find_by_name(Constants.PLASMID)
        self.fragment_type = session.SampleType.find_by_name(Constants.FRAGMENT)
        self.primer_type = session.SampleType.find_by_name(Constants.PRIMER)

    @classmethod
    def new_name(cls):
        return str(uuid4())

    @classmethod
    def _random_record(cls, name: str, cyclic: bool, length_range: Tuple[int, int]):
        length = random.randint(*length_range)
        record = biopython.random_record(
            length, name=name, auto_annotate=True, cyclic=cyclic
        )
        biopython.randomly_annotate(record, feature_length_range=(100, 1000))
        return record

    @classmethod
    def fake_sample(
        cls,
        sample_type: SampleType,
        id: int = None,
        name: str = None,
        properties: dict = None,
    ) -> Sample:
        name = name or cls.new_name()
        sample = sample_type.new_sample(
            name=name, description="", project="", properties=properties
        )
        sample.id = id
        return sample

    def fake_plasmid(
        self,
        id: int = None,
        name: str = None,
        record: biopython.SeqRecord = None,
        length_range: Tuple[int, int] = (4000, 10000),
    ) -> Sample:
        name = name or self.new_name()

        if not record:
            record = self._random_record(
                cyclic=False, name=name, length_range=length_range
            )

        properties = {"Length": len(record.seq)}
        new_sample = self.fake_sample(
            self.plasmid_type, id=id, name=name, properties=properties
        )
        new_sample.record = record
        biopython.make_cyclic([record])
        return new_sample

    def fake_primer(
        self,
        id: int = None,
        name: str = None,
        record: biopython.SeqRecord = None,
        anneal_len_range: Tuple[int, int] = (15, 30),
        overhang_len_range: Tuple[int, int] = (0, 30),
    ) -> Sample:
        name = name or self.new_name()

        if not record:
            anneal_len = random.randint(*anneal_len_range)
            overhang_len = random.randint(*overhang_len_range)
            anneal = biopython.random_record(
                anneal_len, name="anneal", auto_annotate=True, cyclic=False
            )
            overhang = biopython.random_record(
                overhang_len, name="overhang", auto_annotate=True, cyclic=False
            )
            record = overhang + anneal
            record.name = name
        else:
            anneal = record
            overhang = biopython.new_sequence("")
        properties = {
            "Anneal Sequence": str(anneal.seq).upper(),
            "Overhang Sequence": str(overhang.seq).upper(),
            "T Anneal": primer3.calcTm(str(anneal.seq).upper()) - 3,
        }
        new_sample = self.fake_sample(
            self.primer_type, id=id, name=name, properties=properties
        )
        new_sample.record = record
        biopython.make_linear([record])
        return new_sample

    def fake_fragment(
        self,
        id: str = None,
        name: str = None,
        template_record: biopython.SeqRecord = None,
        length_range: Tuple[int, int] = (200, 4000),
    ) -> Sample:
        name = name or self.new_name()
        template = self.fake_plasmid(record=template_record)
        design = primer3plus.Design()
        design.presets.pick_anyway()
        design.presets.template(str(template.record.seq).upper())
        design.presets.product_size(length_range)
        results = design.run_and_optimize(max_iterations=5)
        result = results[0][0]
        size = result["PAIR"]["PRODUCT_SIZE"]
        left = result["LEFT"]["SEQUENCE"]
        right = result["RIGHT"]["SEQUENCE"]
        properties = {
            "Template": template,
            "Forward Primer": self.fake_primer(record=biopython.new_sequence(left)),
            "Reverse Primer": self.fake_primer(record=biopython.new_sequence(right)),
            "Length": size,
        }

        sample = self.fake_sample(
            self.fragment_type, id=id, name=name, properties=properties
        )
        record = template.record[
            result["LEFT"]["location"][0] : result["RIGHT"]["location"][0] + 1
        ]
        assert len(record.seq) == size
        sample.record = record
        biopython.make_linear([record])
        return sample

    def make_fake_library(
        self, n_plasmids: int, n_fragments: int, n_primers: int
    ) -> Tuple[Sample, biopython.SeqRecord]:
        """Generate a list of fake samples (and their SeqRecords)

        :param n_plasmids:
        :param n_fragments:
        :param n_primers:
        :return:
        """
        sample_record_pairs = []
        samples = []
        for _ in range(n_fragments):
            frag = self.fake_fragment()
            template = frag.properties["Template"]
            p1 = frag.properties["Forward Primer"]
            p2 = frag.properties["Reverse Primer"]
            samples += [frag, template, p1, p2]
        for _ in range(n_plasmids):
            samples.append(self.fake_plasmid())
        for _ in range(n_primers):
            samples.append(self.fake_primer())
        for s in samples:
            sample_record_pairs.append((s, s.record))
        return sample_record_pairs
