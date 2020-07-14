import random
import re
from abc import ABC, abstractmethod, abstractproperty
from copy import deepcopy
from functools import partial
from typing import List
from typing import Tuple
from typing import Union, Type
from uuid import uuid4

import validators
from benchlingapi.exceptions import BenchlingAPIException
from benchlingapi.exceptions import InvalidRegistryId
from benchlingapi.models import DNASequence
from Bio.SeqRecord import SeqRecord
from pydent import AqSession
from pydent.models import Sample

from aqbt import bioadapter
from aqbt import biopython
from aqbt import sequence
from aqbt.logger import logger
from aqbt.bioadapter.conversion import seqrecord_to_json

import random

# TODO: move connector and registry to new module? This isn't specific to aquarium


class LabRegistryException(Exception):
    """Generic LabRegistry exceptions."""


class RegistryConnectorABC(ABC):
    @abstractmethod
    def expected_return_type(self) -> Type:
        pass

    @abstractmethod
    def format_registry_id(self, uid: str) -> str:
        pass

    @abstractmethod
    def formatted_registry_id_to_uid(
        self, entity_registry_id: str, uid_pattern: str = r"\d+"
    ):
        pass

    @abstractmethod
    def find_by_url(self, url: str) -> Union[None, DNASequence]:
        pass

    @abstractmethod
    def find(self, uid: str) -> Union[None, DNASequence]:
        pass

    @abstractmethod
    def all(self, limit: int):
        pass


class FakeRegistryConnector(RegistryConnectorABC):
    """Generate a faked DNA registry."""

    def __init__(
        self,
        seq_size_range: Tuple[int, int],
        feature_length_range: Tuple[int, int],
        name: str = "FakeRegistry",
    ):
        """
        Constructor for a connector to a faked DNA registry.

        Usage:

        .. code-block::

            registry = FakeRegistry((1000,5000), (100,2000), name='FakeRegistry')
            registry.all(500)  # generate 500 random sequences
            registry.find("FakeRegistry2")  # get a sequence
            registry.find("FakeRegistry1000", always_return=True)  # generate new seq if it doesn't exist.

        :param seq_size_range:
        :param feature_length_range:
        """
        self.seq_size_range = seq_size_range
        self.feature_length_range = feature_length_range
        self.cache = {}
        self._last_id = 0
        self.name = "FakeRegistry"

    @property
    def expected_return_type(self):
        return dict

    def format_registry_id(self, uid: str) -> str:
        return self.name + str(uid)

    def formatted_registry_id_to_uid(
        self, entity_registry_id: str, uid_pattern: str = r"\d+"
    ):
        pattern = "{name}({uid_pattern})".format(
            name=self.name, uid_pattern=uid_pattern
        )
        m = re.match(pattern, entity_registry_id)
        if m:
            return int(m.group(1))

    def find(self, uid: str, always_return: bool = False) -> Union[None, DNASequence]:
        seq = self.cache.get(uid, None)
        if always_return and seq is None:
            return self._random_seq()
        return seq

    def find_by_url(
        self, url: str, always_return: bool = False
    ) -> Union[None, DNASequence]:
        url_to_seq = {seq["url"] for seq in self.cache.values()}
        seq = url_to_seq.get(url, None)
        if always_return and seq is None:
            return self._random_seq()

    def _random_seq(self):
        length = random.randint(self.seq_size_range[0], self.seq_size_range[1])
        seq = biopython.random_record(
            length, "AGCT", cyclic=[True, False][random.randint(0, 1)]
        )
        biopython.randomly_annotate(seq, feature_length_range=self.feature_length_range)
        self._last_id += 1
        seq = seqrecord_to_json(seq)
        seq["entityRegistryId"] = self.format_registry_id(self._last_id)
        seq["url"] = "www.fakeregistry.org/sequences/{}".format(seq["entityRegistryId"])
        self.cache[seq["entityRegistryId"]] = seq
        return seq

    def all(self, limit: int):
        seqs = []
        for _ in range(limit):
            seqs.append(self._random_seq())
        return seqs


class BenchlingRegistryConnector(RegistryConnectorABC):
    def __init__(self, api, initials, schema, prefix, folder_id, registry_id):
        self.api = api
        self.organization_initials = initials
        self.schema = schema
        self.prefix = prefix
        self.folder_id = folder_id
        self.folder = self.api.Folder.find(self.folder_id)
        self.registry = self.api.Registry.find(registry_id)
        self.convert = partial(
            bioadapter.convert,
            benchling_folder_id=self.folder_id,
            benchling_session=self.api,
        )
        self.logger = logger(self)

    @property
    def expected_return_type(self):
        return DNASequence

    def format_registry_id(self, uid: str) -> str:
        """Format a unique id to a Benchling registry id.

        :param uid: unique id
        :return: Bencling registry id
        """
        return "{prefix}{init}{uid}".format(
            prefix=self.prefix, uid=uid, init=self.organization_initials
        )

    def formatted_registry_id_to_uid(
        self, entity_registry_id: str, uid_pattern: str = r"\d+"
    ) -> str:
        pattern = "{prefix}{init}({uid_pattern})".format(
            prefix=self.prefix, init=self.organization_initials, uid_pattern=uid_pattern
        )
        m = re.match(pattern, entity_registry_id)
        if m:
            return int(m.group(1))

    def find_by_url(self, url: str) -> Union[None, DNASequence]:
        """Find a DNASequence by its url or sharelink.

        :param url: url or share link
        :return: DNASequence or None
        """
        dna = None
        if validators.url(url) is True:
            try:
                dna = self.api.DNASequence.from_share_link(url)
            except BenchlingAPIException:
                dna = None
            if not dna:
                self.logger.debug("Could not find DNA using '{}'".format(url))
        return dna

    def find(self, uid: str) -> Union[None, DNASequence]:
        """Find a DNASequence in the Benchling registry by its unique id, which
        will get converted to a Benchling registry id.

        :param uid: unique id
        :return: DNASequence or None
        """
        rid = self.format_registry_id(uid)
        try:
            return self.api.DNASequence.find_in_registry(rid)
        except InvalidRegistryId:
            self.logger.debug("Could not find '{}' in registry".format(rid))
            return None

    # def might_by_registered(self, seq: DNASequence):
    #     if seq.entity_registry_id:
    #         return True
    #     return self.is_registered(seq)

    def register_sequence(
        self, seq: DNASequence, uid: str, overwrite: bool = True, do_raise=False
    ):
        """Register a DNASequence to Benchling.

        :param seq: the DNASequence to register
        :param uid: the unique id (will be converted to benchling registry id)
        :param overwrite: if True (default) and if the entity exists in the registry,
            the existing entity will be updated from the provided sequence. Else,
            the update will be ignored if do_raise == False, or an `InvalidRegistryId`
            will be raised if do_raise == True
        :param do_raise: see `overwrite` argument
        :return: the registered DNASequence
        """
        entity_registry_id = self.format_registry_id(uid)
        existing = self.find(uid)
        if existing:
            if overwrite:
                self.logger.info("Updating existing sequence {}".format(existing.id))
                self.api.DNASequence.update_model(existing.id, seq.update_json())
            elif do_raise:
                self.logger.error("Cannot overwrite.")
                raise InvalidRegistryId("Sequence already registered.")
            else:
                self.logger.info(
                    "Registration ignored, returning existing sequence {}".format(
                        existing.entity_registry_id
                    )
                )
            seq = existing
        else:
            if not hasattr(seq, "id") or seq.id is None:
                seq.folder_id = self.folder_id
                seq.save()
            seq.set_schema(self.schema)
            seq.register_with_custom_id(entity_registry_id)
            self.logger.info("Registered {}".format(entity_registry_id))
        return seq

    def all(self, limit: int):
        return self.api.DNASequence.all(registry_id=self.registry.id, limit=limit)


def make_fake_seq_record(
    mn_bases: int,
    mx_bases: int,
    cyclic_choices: Tuple[bool] = (True, False),
    feature_size_range: Tuple = (100, 1500),
) -> SeqRecord:
    n_bases = random.randint(mn_bases, mx_bases)
    cyclic = random.sample(cyclic_choices, k=1)[0]
    rec = biopython.random_record(n_bases, name=None, cyclic=cyclic, auto_annotate=True)
    biopython.randomly_annotate(rec, feature_length_range=feature_size_range)
    biopython.randomly_annotate(rec, feature_length_range=feature_size_range)
    return rec


class RetrievalPriorities:

    PRIMER = "primer"
    CACHE = "cache"
    REGISTRY = "benchling_registry"
    URL = "url_or_sharelink"
    UNANNOTATED_PROPERTY = "unannotated_property"
    VALID_PROPERTIES = [PRIMER, CACHE, REGISTRY, URL, UNANNOTATED_PROPERTY]


# TODO: move config to init
CONFIG = {
    "primer_sample_type": "Primer",
    "fragment_sample_type": "Fragment",
    "plasmid_sample_type": "Plasmid",
    "default_sequence_field_names": ("Sequence",),
    "primer_sequence_field_names": ("Overhang Sequence", "Anneal Sequence"),
    "fragment_template_field_name": "Template",
    "fragment_primer_field_names": ("Forward Primer", "Reverse Primer"),
}


class LabDNARegistry:
    """Establishes connection between Aquarium and a DNA Registry."""

    def __init__(self, connector: BenchlingRegistryConnector, aqsession: AqSession):
        self.session = aqsession

        self._primer_type = self.session.SampleType.find_by_name(
            CONFIG["primer_sample_type"]
        )
        self._fragment_type = self.session.SampleType.find_by_name(
            CONFIG["fragment_sample_type"]
        )
        self._plasmid_type = self.session.SampleType.find_by_name(
            CONFIG["plasmid_sample_type"]
        )
        assert self._primer_type
        assert self._fragment_type
        assert self._plasmid_type
        self.connector = connector
        self._using_cache = False
        self._registry_cache = {}
        self.logger = logger(self)

    # TODO: this is specific for benchling
    @property
    def benchling(self):
        return self.connector.api

    @property
    def using_cache(self):
        return self._using_cache

    @using_cache.setter
    def using_cache(self, b):
        if b and not self._using_cache:
            self.refresh_cache()
        self._using_cache = b

    def use_cache(self, lim: int = None):
        if not self._using_cache:
            self.refresh_cache(lim)
        self._using_cache = True

    def refresh_cache(self, limit: int = None):
        self._registry_cache = {}
        for dna in self.connector.all(limit=limit):
            if dna.entity_registry_id:
                self.add_to_cache(dna)

    # TODO: this is specific for benchling
    def _make_fake_benchling_sequence(
        self,
        mn_bases: int = 500,
        mx_bases: int = 10000,
        cyclic_choices: Tuple[bool] = (True, False),
        feature_size_range: Tuple = (100, 1500),
        entity_registry_id: str = None,
    ) -> DNASequence:
        dna = bioadapter.convert(
            make_fake_seq_record(
                mn_bases, mx_bases, cyclic_choices, feature_size_range
            ),
            to="DNASequence",
            benchling_session=self.benchling,
            benchling_folder_id="not_a_folder",
        )
        if not entity_registry_id:
            entity_registry_id = str(uuid4())
        dna.entity_registry_id = entity_registry_id
        return dna

    def use_fake_cache(self, n_seqs: int = None, entity_registry_ids: List[str] = None):
        self._registry_cache = {}
        if entity_registry_ids:
            for eid in entity_registry_ids:
                self.add_to_cache(
                    self._make_fake_benchling_sequence(entity_registry_id=eid)
                )
        else:
            for i in range(n_seqs):
                self.add_to_cache(self._make_fake_benchling_sequence())
        self._using_cache = True

    def add_to_cache(self, dna: DNASequence):
        if dna.entity_registry_id:
            self._registry_cache[dna.entity_registry_id] = dna
        else:
            raise ValueError("DNASequence is missing entity_registry_id")

    def find_in_cache(self, sample: Sample) -> Union[None, DNASequence]:
        if not sample.id:
            raise ValueError("Sample has not id")
        uid = self.connector.format_registry_id(sample.id)
        return self._registry_cache.get(uid, None)

    def find_in_registry(self, sample: Sample) -> Union[None, DNASequence]:
        if not sample.id:
            raise ValueError("Sample has not id")
        return self.connector.find(sample.id)

    @staticmethod
    def _safe_get_sample_property(sample: Sample, keys: Tuple[str]):
        properties = sample.properties
        property = None
        for key in keys:
            try:
                property = properties.get(key, None)
                break
            except Exception as e:
                pass
        return property

    def find_by_url(
        self, sample: Sample, keys: Tuple[str] = None
    ) -> Union[None, DNASequence]:
        """Find a DNASequence from an Aquarium Sample a url located int its
        property keys. By default 'Sequence' is the only key used.

        :param sample: the Aquarium sample
        :param keys: list of keys to find the url
        :return: DNASequence or None
        """
        if keys is None:
            keys = CONFIG["default_sequence_field_names"]
        url = self._safe_get_sample_property(sample, keys)
        return self.connector.find_by_url(url)

    def find_by_property(
        self, sample: Sample, keys: Tuple[str] = None
    ) -> Union[None, DNASequence]:
        """Find a DNASequence from an Aquarium Sample its sequence located int
        its property keys. By default 'Sequence' is the only key used.

        :param sample: the Aquarium sample
        :param keys: list of keys to find the sequence str
        :return: DNASequence or None
        """
        if keys is None:
            keys = (CONFIG["default_sequence_field_names"],)
        sequence_str = self._safe_get_sample_property(sample, keys)
        if sequence.dna_like(sequence_str):
            return biopython.new_sequence(
                sequence_str, name=sample.name, auto_annotate=True
            )

    def get_primer_sequence(
        self, sample: Sample, keys: Tuple[Union[str, Tuple[str, str]]] = None,
    ) -> SeqRecord:
        """Return the primer sequence as a SeqRecord sequence from an Aquarium
        sample using its properties.

        :param sample: the Aquarium sample
        :param keys: ordered tuple of keys to retrieve the sequences from. SeqRecords
            will be created from each part and concatenated. If provided with a tuple
            of tuples, the parts will automatically be annotated with a feature
            from the second entry ((key1, name1), (key2, name2), ...).
        :return: the SeqRecord
        """
        if keys is None:
            keys = CONFIG["primer_sequence_field_names"]
        if sample.sample_type_id == self._primer_type.id:
            record = None
            for key in keys:
                name = key
                if isinstance(key, tuple):
                    key, name = key
                elif not isinstance(key, str):
                    raise ValueError(
                        "Keys must contain either strings or (str, str) tuples."
                    )
                part = sample.properties.get(key, "") or ""
                record_part = biopython.new_sequence(part, name)
                if record is None:
                    record = record_part
                else:
                    record += record_part
            record.name = sample.name
            record.id = sample.id
            biopython.annotate(record, sample.name)
            return record
        raise TypeError(
            "Sample {} is not a PrimerType. {}".format(
                sample.name, sample.sample_type.name
            )
        )

    def get_sequence(
        self,
        sample: Sample,
        return_type="DNASequence",
        ignore_registry: bool = False,
        priority: List[str] = None,
    ):
        """Return the Aquarium sequence as a DNASequence (or specified type).

        1. look for 'seq' attribute
        2. look sequence from the Benchling Registry using its entity registry id
        3. look for str in sample.properties['Sequence']. Try to convert the string
            into a SeqRecord. Else, use the string as a weburl for finding the sequence.
        4. Convert the sequence to the specified type.

        :param sample: Aquarium Sample to get the sequence
        :param return_type: type to convert (default Benchling style DNASequence)
        :param ignore_registry: if True, will ignore retrieving from the registry
        :param priority:
        :return:
        """

        if priority is None:
            priority = [
                RetrievalPriorities.CACHE,
                RetrievalPriorities.REGISTRY,
                RetrievalPriorities.URL,
                RetrievalPriorities.UNANNOTATED_PROPERTY,
            ]

        fxns = {
            RetrievalPriorities.CACHE: self.find_in_cache,
            RetrievalPriorities.PRIMER: self.get_primer_sequence,
            RetrievalPriorities.REGISTRY: self.find_in_registry,
            RetrievalPriorities.UNANNOTATED_PROPERTY: self.find_by_property,
            RetrievalPriorities.URL: self.find_by_url,
        }

        for p in priority:
            fxn = fxns[p]
            seq = fxn(sample)
            if seq is not None:
                break

        if not seq:
            self.logger.error("Could not find sequence.")
            return None

        seq = self.connector.convert(seq, to=return_type)
        return seq
        #
        # self.logger.debug("Attempting to find sequence.")
        # seq = None
        # if hasattr(sample, "seq"):
        #     seq = sample.seq
        #     msg = "Found sequence on instance 'seq' attribute."
        #
        # if not ignore_registry:
        #     if self.use_cache and not seq and sample.id:
        #         seq = self.find_in_cache(sample)
        #         msg = "Found in registry cache."
        #     if not seq and sample.id:
        #         seq = self.connector.find(sample.id)
        #         msg = "Found on registry."
        # if not seq:
        #     try:
        #         seq = sample.properties.get("Sequence", None)
        #     except Exception:
        #         seq = None
        #
        # if isinstance(seq, str):
        #     if sequence.dna_like(seq):
        #         seq = biopython.new_sequence(seq, name=sample.name, auto_annotate=True)
        #         msg = "Created new SeqRecord from 'Sequence' sample properties."
        #     else:
        #         try:
        #             seq = self.connector.api.DNASequence.from_share_link(seq)
        #         except Exception:
        #             seq = None
        #         msg = "Found sequence using benchling sharelink."
        # if not seq:
        #     self.logger.error("Could not find sequence.")
        #     return None
        # self.logger.debug(msg)
        # seq = self.connector.convert(seq, to=return_type)
        # return seq

    def is_registered(self, sample: Sample) -> bool:
        """Check that an Aquarium sample has a registered DNA sequence."""
        return self.get_sequence(sample) is not None

    def fast_is_registered(self, sample: Sample) -> bool:
        """Quickly check if an Aquarium sample has a registered DNA sequence
        by checking the cached sequences."""
        if not self._using_cache:
            raise ValueError(
                "`registry.use_cache` must be True to use '{}'."
                " Set it using `registry.use_cache = True`".format(
                    self.fast_is_registered.__name__
                )
            )
        found = self.find_in_cache(sample)
        if found:
            return found
        return False

    def register(
        self,
        sample: Sample,
        seq: Union[DNASequence, SeqRecord],
        overwrite=True,
        do_raise=True,
    ):
        """Register an Aquarium sequence sample."""
        if not sample.id:
            raise LabRegistryException(
                "Cannot registry sample because" " it has not been saved in Aquarium."
            )
        if seq is None:
            seq = self.get_sequence(sample, return_type="DNASequence")
        if not seq:
            raise LabRegistryException(
                "Cannot registry sample because it has no sequence attached."
            )
        if isinstance(seq, SeqRecord):
            seq = self.connector.convert(seq, to="DNASequence")
        if hasattr(seq, "primers") and seq.primers:
            seq.primers = []
        registered = self.connector.register_sequence(
            seq, sample.id, overwrite=overwrite, do_raise=do_raise
        )
        if registered and self.using_cache:
            self.add_to_cache(registered)
        return registered

    def copy(self) -> "LabDNARegistry":
        copied = self.__class__(connector=self.connector, aqsession=self.session())
        copied._registry_cache = deepcopy(self._registry_cache)
        copied._using_cache = self._using_cache
        return copied

    def pcr_products(self, fragment: Sample) -> List[SeqRecord]:
        primers = [
            self.get_primer_sequence(fvname)
            for fvname in CONFIG["fragment_primer_field_names"]
        ]
        template = self.get_sequence(
            fragment.properties[CONFIG["fragment_template_field_name"]]
        )
        template_record = self.connector.convert(template, to="DNASequence")

        product_records = biopython.pcr_amplify(
            primers, template_record, cyclic=template.is_circular, name=fragment.name
        )
        return product_records
