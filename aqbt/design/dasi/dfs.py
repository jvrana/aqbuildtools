import hashlib
from typing import List
from typing import Union
from uuid import uuid4

import inflection
import pandas as pd
from benchlingapi.models import DNASequence
from Bio.SeqRecord import SeqRecord
from pydent import models
from tqdm import tqdm

from aqbt import bioadapter
from aqbt import biopython
from aqbt.aquarium.faker import FakeSampleGenerator
from aqbt.aquarium.pydent_utils import Constants as C
from aqbt.aquarium.registry import LabDNARegistry
from aqbt.utils import chunkify


def seq_sha1(seq: str) -> str:
    """Convert sequence string into a hash."""
    return hashlib.sha1(seq.strip().upper().encode()).hexdigest()


# TODO: initialize with config
CONFIG = {
    "plasmid_sample_type": "Plasmid",
    "fragment_sample_type": "Fragment",
    "primer_sample_type": "Primer",
    "fragment_object_types": ["Fragment Stock"],
    "plasmid_object_types": ["Plasmid Glycerol Stock"],
    "linear_types": ["Primer", "Fragment"],
    "cyclic_types": ["Plasmid"],
}


# TODO: use Aquarium config to determine connections
class KlavinsLabDnaDb:
    """Generates the inventory database for DASi."""

    DASI_DEFAULT_AQ_TIMEOUT = 120  #: timout (s) for Aquarium API
    LIMS_ID = "LIMS_ID"  #: lims id key
    MAX_AQ_ID = 10 ** 7  #: maximum expected Aquarium ID
    BROWSER_MAX_CHUNK_SIZE = 100  #: maximum size the browser can use to cache things
    EXPECTED_COLUMNS = (
        "sample_id",
        "record",
        "sample_type",
        "is_available",
        "sequence_hash",
        "record_uuid",
    )  #: required columns
    UUID = "UUID"
    ALL_COLUMNS = [
        "sample_id",  # the sample id of the Aquarium sample
        "benchling_sequence",  # the Benchling sequence (benchling.DNASequence)
        "record",  # the biopython sequence (Bio.SeqRecord)
        "entity_registry_id",  # the entity_registry_id
        "is_circular",  # whether the topology of the sequence
        "is_available",  # whether there is available inventory for the sample
        "sample",  # the Aquarium Sample instance
        "sample_type",  # the name of the SampleType
        "sequence_hash",  # a SHA1 hash of the sequence for sequence comparison
        "record_uuid",  # the unique ID for the record / sample
    ]

    def __init__(self, registry: LabDNARegistry):
        self.registry = registry

        self.plasmid_type = self.session.SampleType.find_by_name(
            CONFIG["plasmid_sample_type"]
        )
        self.fragment_type = self.session.SampleType.find_by_name(
            CONFIG["fragment_sample_type"]
        )
        self.primer_type = self.session.SampleType.find_by_name(
            CONFIG["primer_sample_type"]
        )

        assert self.plasmid_type
        assert self.fragment_type
        assert self.primer_type

        self.valid_fragment_object_ids = [
            self.session.ObjectType.find_by_name(n).id
            for n in CONFIG["fragment_object_types"]
        ]
        self.valid_plasmid_object_ids = [
            self.session.ObjectType.find_by_name(n).id
            for n in CONFIG["plasmid_object_types"]
        ]

        self.valid_linear_types = CONFIG["linear_types"]
        self.valid_cyclic_types = CONFIG["cyclic_types"]

        self.df = None

    @property
    def session(self):
        return self.registry.session

    @staticmethod
    def _series(**kwargs):
        keys = list(kwargs.keys())
        values = list(kwargs.values())
        return pd.Series(values, index=keys)

    def new_row(
        self,
        sample: models.Sample,
        record: SeqRecord,
        is_available: bool,
        is_circular: bool = None,
        sample_id: int = None,
        benchling_sequence: DNASequence = None,
        entity_registry_id: str = None,
        sample_type: str = None,
        sequence_hash: str = None,
        record_uuid: str = None,
    ) -> pd.Series:
        if sample_id is None:
            sample_id = sample.id

        if record_uuid is None:
            record_uuid = str(uuid4())

        if sequence_hash is None:
            self._make_seq_hash(record)

        if entity_registry_id is None:
            self.registry.connector.format_registry_id(sample_id)

        if is_circular is None:
            is_circular = biopython.is_circular(record)

        if sample_type is None:
            sample_type = sample.sample_type.name
        return self._series(
            sample_id=sample_id,
            benchling_sequence=benchling_sequence,
            record=record,
            entity_registry_id=entity_registry_id,
            is_circular=is_circular,
            is_available=is_available,
            sample=sample,
            sample_type=sample_type,
            sequence_hash=sequence_hash,
            record_uuid=record_uuid,
        )

    def build_fake_inventory_df(
        self, n_plasmid: int, n_fragment: int, n_primers: int
    ) -> pd.DataFrame:
        """Build a faked inventory df for testing / debugging purposes.

        :param registry: klavins lab registry
        :param n_plasmid: number of plasmids to generate
            (does not include primers from generated from fake_fragments)
        :param n_fragment: number of fragments to generate, including new templates and primers
        :param n_primers: number of primers to generate
            (does not include primers from generated from fake_fragments)
        :return: fake inventory df
        """
        gen = FakeSampleGenerator(self.registry.session())
        sample_and_records = gen.make_fake_library(n_plasmid, n_fragment, n_primers)

        series = []
        for i, (sample, record) in enumerate(sample_and_records):
            sample.id = i
            series.append(self.new_row(sample=sample, record=record, is_available=True))
        df = pd.DataFrame(series)
        self.post_process_df(df)
        self.validate_df(df)
        return df

    @staticmethod
    def sanitize_record(rec: SeqRecord) -> SeqRecord:
        rec.name = inflection.parameterize(rec.name)
        return rec

    @staticmethod
    def _set_row_topology(row: pd.Series, is_circular: bool):
        row["is_circular"] = is_circular
        if is_circular:
            biopython.make_cyclic([row["record"]])
            row["benchling_sequence"].is_circular = True
        else:
            biopython.make_linear([row["record"]])
            row["benchling_sequence"].is_circular = False

    def validate_df(
        self,
        df: pd.DataFrame,
        ignore: bool = False,
        auto_correct_topologies: bool = True,
    ):
        """

        :param df:
        :param ignore:
        :param auto_correct_topologies: if True, will autocorrect the SeqRecord and Benchling DNASequence
            topologies for the dataframe.
        :return:
        """
        # check for expected columns
        for c in self.EXPECTED_COLUMNS:
            if c not in df.columns:
                raise ValueError("Column '{}' is missing from DF.".format(c))

        # check topologies
        rows = []
        for i, row in df.iterrows():
            if not row["is_circular"] == biopython.is_circular(row["record"]):
                raise ValueError("Record and 'is_circular' key does not match.")
            if row["sample_type"] in self.valid_linear_types:
                if row["is_circular"] is True:
                    msg = "{} (sample_id={}, row={}) cannot be circular".format(
                        row["sample_type"], row["sample_id"], i
                    )
                    if auto_correct_topologies:
                        print("Correcting topology. ", end="")
                        print(msg)
                        self._set_row_topology(row, is_circular=False)
                        rows.append(row)
                    elif ignore:
                        print(msg)
                    else:
                        raise ValueError(msg)
                else:
                    rows.append(row)
            elif row["sample_type"] in self.valid_cyclic_types:
                if row["is_circular"] is False:
                    msg = "{} (sample_id={}, row={}) cannot be linear".format(
                        row["sample_type"], row["sample_id"], i
                    )
                    if auto_correct_topologies:
                        print("Correcting topology. ", end="")
                        print(msg)
                        self._set_row_topology(row, is_circular=True)
                        rows.append(row)
                    elif ignore:
                        print(msg)
                    else:
                        raise ValueError(msg)
                else:
                    rows.append(row)
            rows.append(row)
        return pd.DataFrame(rows)

    @classmethod
    def annotate_record_with_lims_id(cls, record: SeqRecord, lims_id: Union[str, int]):
        """Adds a LIMS_ID to the record annotations dictionary."""
        record.annotations[cls.LIMS_ID] = lims_id

    def post_process_df(self, df: pd.DataFrame):
        """Post processing for inventory dataframes.

        :param df:
        :return:
        """
        # drop duplicates
        df.drop_duplicates(inplace=True, subset=["sample_id"])

        # filter out invalid sample types
        df = df.loc[
            lambda df: df["sample_type"].isin(
                self.valid_cyclic_types + self.valid_linear_types
            )
        ]

        for _, row in df.iterrows():
            record = row["record"]

            # add information to SeqRecord.annotations
            self.annotate_record_with_lims_id(record, row["sample_id"])
            record.annotations["uuid"] = row["record_uuid"]
            record.name = row["sample"].name
            # sanitize record name
            self.sanitize_record(record)

            # annotate with topology
            if row["is_circular"]:
                biopython.make_cyclic(record)
            else:
                biopython.make_linear(record)

            # check topologies
            if row["is_circular"]:
                if not row["sample_type"] == C.PLASMID:
                    raise ValueError(
                        "{} must be {}".format(row["sample_type"], C.PLASMID)
                    )
            else:
                if not row["sample_type"].strip() in [C.PRIMER, C.FRAGMENT]:
                    raise ValueError(
                        "{} must be {}".format(
                            row["sample_type"], [C.PRIMER, C.FRAGMENT]
                        )
                    )

            # attach record to sample
            row["sample"].record = record

    def build(
        self, dna_limit: int = None, primer_limit: int = None, ignore: bool = False
    ) -> pd.DataFrame:
        primer_df = self.build_primer_df(limit=primer_limit)
        primer_df = self.validate_df(primer_df, ignore=ignore)

        plasmid_and_fragment_df = self.build_inventory_df(limit=dna_limit)
        plasmid_and_fragment_df = self.validate_df(
            plasmid_and_fragment_df, ignore=ignore
        )

        df = primer_df.append(plasmid_and_fragment_df)
        df = self.validate_df(df, ignore=ignore)

        self.post_process_df(df)
        self.df = df
        return df

    def build_primer_df(self, limit: int = None) -> pd.DataFrame:
        """Build primer df."""

        list_of_series = []
        columns = self.ALL_COLUMNS[:]
        with self.session.with_cache(
            timeout=self.DASI_DEFAULT_AQ_TIMEOUT, using_models=True
        ) as sess:
            primer_type = sess.SampleType.find_by_name(C.PRIMER)
            print("collecting all primer samples from Aquarium")

            query = {"sample_type_id": primer_type.id}
            if limit:
                primers = sess.Sample.last(limit, query=query)
            else:
                primers = sess.Sample.where(query)

            print("collecting all primer properties from Aquarium")
            primer_records = []
            for _primers in chunkify(primers, 1000):
                sess.browser.get(
                    _primers,
                    {
                        "field_values": {"field_type": "allowable_field_types"},
                        "sample_type": "field_types",
                    },
                )
                for primer in tqdm(_primers, desc="created sequences for primers"):
                    record = None
                    try:
                        record = self.registry.get_primer_sequence(primer)
                    except Exception as e:
                        print(
                            "Failed to retrieve sequence for {} because {}".format(
                                primer.name, str(e)
                            )
                        )
                    if record and str(record.seq):
                        record.id = "{}__{}".format(primer.id, primer.name)
                        primer_records.append(record)
                        self.annotate_record_with_lims_id(record, primer.id)
                        list_of_series.append(
                            pd.Series(
                                [
                                    primer.id,
                                    None,
                                    record,
                                    None,
                                    False,
                                    True,
                                    primer,
                                    primer_type.name,
                                    seq_sha1(str(record.seq)),
                                    str(uuid4()),
                                ],
                                index=columns,
                            )
                        )
        return pd.DataFrame(list_of_series, columns=columns)

    def build_inventory_df(self, limit: int = None) -> pd.DataFrame:
        """Construct a pd.DataFrame with the following columns.

        sample_id | benchling_sequence | record | entity_registry_id | is_circular |
        is_available | sample | sample_type

        :param registry:
        :param limit:
        :return:
        """
        self.registry.use_cache(limit)

        dnas = list(self.registry._registry_cache.values())
        df1 = self._build_record_df(dnas)
        sample_ids = list(df1["sample_id"])
        if limit:
            sample_ids = sample_ids[-limit:]
        df2 = self._build_sample_df(sample_ids)

        df3 = df1.set_index("sample_id").join(df2.set_index("sample_id"))
        df3["sample_id"] = df3.index
        return df3

    @staticmethod
    def _make_seq_hash(record: SeqRecord) -> str:
        return seq_sha1(str(record.seq))

    def _build_record_df(self, dnas: List[DNASequence]) -> pd.DataFrame:
        """Construct a DataFrame of all dna sequences from Benchling.

        :param registry:
        :param dnas:
        :return:
        """
        # dna df
        list_of_series = []
        for dna in dnas:
            sid = self.registry.connector.formatted_registry_id_to_uid(
                dna.entity_registry_id
            )
            record = bioadapter.convert(dna, to="SeqRecord")
            if str(record.seq):
                list_of_series.append(
                    pd.Series(
                        [
                            sid,
                            dna,
                            record,
                            dna.entity_registry_id,
                            dna.is_circular,
                            self._make_seq_hash(record),
                            str(uuid4()),
                        ],
                        index=[
                            "sample_id",
                            "benchling_sequence",
                            "record",
                            "entity_registry_id",
                            "is_circular",
                            "sequence_hash",
                            "record_uuid",
                        ],
                    )
                )
        return pd.DataFrame(list_of_series)

    def fragment_is_available(
        self, sample: models.Sample
    ) -> Union[None, List[models.Item]]:
        ok_items = []
        for i in sample.items:
            if (
                i.location != C.DELETED
                and i.object_type_id in self.valid_fragment_object_ids
            ):
                ok_items.append(i)
        if ok_items:
            return True
        return False

    def plasmid_is_available(
        self, sample: models.Sample
    ) -> Union[None, List[models.Item]]:
        ok_items = []
        for i in sample.items:
            if (
                i.location != C.DELETED
                and i.object_type_id in self.valid_plasmid_object_ids
            ):
                ok_items.append(i)
        return ok_items

    def dna_is_available(self, sample: models.Sample):
        """Checks if there the sample has 'available items'."""
        if sample is None:
            raise ValueError()
        if sample.sample_type_id == self.plasmid_type.id:
            return self.plasmid_is_available(sample)
        elif sample.sample_type_id == self.fragment_type.id:
            return self.fragment_is_available(sample)

    @classmethod
    def _safe_sample_ids(cls, sample_ids: List[int]) -> List[int]:
        return [sid for sid in sample_ids if sid and 0 < sid < cls.MAX_AQ_ID]

    def _build_sample_df(self, sample_ids: List[int]):
        """Construct the 'Sample' DataFrame with columns `["sample_id",
        "sample", "sample_type", "is_available"]`

        DataFrame Keys:

        ::
        sample_id: the Aquarium sample_id
        sample_type: name of the Aquarium sample_type
        sample: the Aquarium Sample instance
        is_available: whether the Sample has a valid inventory type

        :param session:
        :param sample_ids: list of Aquarium sample_ids
        :return:
        """
        sample_ids = self._safe_sample_ids(sample_ids)

        list_of_series = []

        with self.session.with_cache(
            timeout=self.DASI_DEFAULT_AQ_TIMEOUT, using_models=True
        ) as sess:
            print("retrieving samples...", end=" ")

            registered_samples = sess.Sample.where({"id": sample_ids})
            sess.browser.get(registered_samples, "sample_type")
            print(len(registered_samples))
            print("retrieving items...", end=" ")

            for _samples in chunkify(registered_samples, self.BROWSER_MAX_CHUNK_SIZE):
                sess.browser.get(_samples, "items")
            for sample in registered_samples:
                is_available = False
                available_items = self.dna_is_available(sample)
                if available_items:
                    is_available = True
                list_of_series.append(
                    pd.Series(
                        [sample.id, sample, sample.sample_type.name, is_available],
                        index=["sample_id", "sample", "sample_type", "is_available"],
                    )
                )
        return pd.DataFrame(list_of_series)

    @property
    def available_df(self):
        return self.df[self.df["is_available"] is True]
