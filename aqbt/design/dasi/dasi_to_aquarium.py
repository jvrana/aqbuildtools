import pickle
from os.path import abspath
from os.path import dirname
from os.path import isfile
from os.path import join

from dasi import LibraryDesign

from aqbt.bioadapter import convert
from aqbt.design.dasi.dfs import KlavinsLabDnaDb


here = abspath(dirname("."))
DEFAULT_FILEPATH = join(here, "inventory.df")


def inventory_df(registry=None, force_new: bool = False, filepath: str = None):
    if filepath is None:
        filepath = DEFAULT_FILEPATH
    df = None
    if isfile(filepath) and not force_new:
        with open(filepath, "rb") as f:
            try:
                return pickle.load(f)
            except Exception as e:
                print("Could not load pickled file.")
                print(e)
    if df is None:
        db = KlavinsLabDnaDb(registry)
        db.build(dna_limit=None, primer_limit=None)

        df = db.df[db.df["is_available"] == True]
        df.drop("benchling_sequence", axis=1, inplace=True)
        with open(filepath, "wb") as f:
            pickle.dump(df, f)

    # clean
    df["trashed"] = [s.name.lower().startswith("trashed") for s in df["sample"]]
    df = df[df["is_available"] == True]
    df = df[~df["trashed"]]
    return df


def design_from_benchling_links(links, registry, df):
    seqs = [registry.connector.api.DNASequence.from_share_link(link) for link in links]

    records = [convert(seq, to="SeqRecord") for seq in seqs]

    df = df[df["is_available"] is True]
    df = df[~df["trashed"]]

    def _get_records(st):
        return list(df[df["sample_type"] == st]["record"])

    fragments = _get_records("Fragment")
    plasmids = _get_records("Plasmid")
    primers = _get_records("Primer")

    design = LibraryDesign()
    design.FAVOR_SHARED_SEQUENCES = 2.0
    design.add_materials(
        primers=primers, fragments=fragments, templates=plasmids, queries=records
    )
    design.logger.set_level("DEBUG")
    design.run(n_paths=1)

    return design
