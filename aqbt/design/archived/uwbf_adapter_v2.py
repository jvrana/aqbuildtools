import pandas as pd
from lobio.klavinslab.pydent_utils import Constants as C
from lobio.klavinslab.registry import KlavinsLabRegistry

from .dfs import build_df


class DesignFactory:

    # TODO: add debug_n_fragments
    def __init__(
        self,
        registry: KlavinsLabRegistry,
        directory: str = ".",
        debug_mode: bool = False,
        debug_n_seqs: int = 100,
        debug_primer_lim: int = 100,
        inventory_df: pd.DataFrame = None,
        ignore_failed_sequences: bool = False,
    ):
        """

        :param registry:
        :param directory:
        :param debug_mode:
        :param debug_n_seqs:
        :param debug_primer_lim:
        :param inventory_df: optional inventory df. If None, will produce a new inventory_df.
        """

        # settings
        self.dir = directory  #: the directory to save
        self.do_load = True  #: if True, will load SpanCost bytes, if available
        self.do_save = True  #: if True, will save SpanCost bytes, if available
        self.debug_n_seqs = (
            debug_n_seqs  #: number of sequences to run while in debug mode
        )
        self.debug_mode = debug_mode  #: run the design in debug mode, generating random Samples, Items, and sequences

        self.registry = registry

        # make inventory df
        if inventory_df is not None:
            self._inventory_df = inventory_df
        elif self.debug_mode:
            self._inventory_df = build_df(
                self.registry,
                dna_limit=debug_n_seqs,
                primer_limit=debug_primer_lim,
                ignore=ignore_failed_sequences,
            )
        else:
            self._inventory_df = build_df(
                self.registry,
                dna_limit=None,
                primer_limit=None,
                ignore=ignore_failed_sequences,
            )

        # default sample types
        self.fragment_type = self.registry.session.SampleType.find_by_name(C.FRAGMENT)
        self.plasmid_type = self.registry.session.SampleType.find_by_name(C.PLASMID)
        self.primer_type = self.registry.session.SampleType.find_by_name(C.PRIMER)

    @property
    def inventory_df(self) -> pd.DataFrame:
        """
        A pandas.DataFrame that contains inventory information
        for all the samples that are available. the following columns:

        ::

            ALL_COLUMNS = [
                "sample_id",            # the sample id of the Aquarium sample
                "benchling_sequence",   # the Benchling sequence (benchling.DNASequence)
                "record",               # the biopython sequence (Bio.SeqRecord)
                "entity_registry_id",   # the entity_registry_id
                "is_circular",          # whether the topology of the sequence
                "is_available",         # whether there is available inventory for the sample
                "sample",               # the Aquarium Sample instance
                "sample_type",          # the name of the SampleType
                "sequence_hash",        # a SHA1 hash of the sequence for sequence comparison
                "record_uuid",          # the unique ID for the record / sample
            ]

        :return: the dataframe
        """
        return self._inventory_df[self._inventory_df["is_available"] == True]
