from json.decoder import JSONDecodeError
from typing import Union

from tqdm import tqdm
from tqdm import tqdm_gui
from tqdm import tqdm_notebook

from aqbt.aquarium.parsers.pgrr_plasmid import parse_name
from aqbt.bioadapter import convert


def register_all(
    registry, tqdm_class: Union[tqdm_notebook, tqdm_gui, tqdm] = tqdm_notebook
):
    benchling = registry.benchling
    session = registry.session
    registry.using_cache = True
    with session.with_cache(timeout=60) as sess:
        plasmid_type = sess.SampleType.find_by_name("Plasmid")
        plasmids = sess.Sample.where({"sample_type_id": plasmid_type.id})

        iterator = tqdm_class(plasmids, desc="collecting pGRR plasmids")
        unregistered = []
        for p in iterator:
            try:
                is_registered = registry.fast_is_registered(p)
            except JSONDecodeError:
                is_registered = False
            if not is_registered:
                seq = parse_name(p.name)
                if seq:
                    unregistered.append((p, seq))

        iterator = tqdm_class(unregistered, desc="registering pGRR plasmids")
        for p, seq in iterator:
            seq = parse_name(p.name)
            if seq:
                dna = convert(
                    seq,
                    to="DNASequence",
                    benchling_session=benchling,
                    benchling_folder_id=registry.connector.folder.id,
                )
                dna.primers = []
                print('Registering "{}"'.format(dna.name))
                dna.merge(on=["name", "folder_id"])
                registry.register(sample=p, seq=dna)

                print("Registration successful!")
                print(dna.web_url)


# def test_register_all():
#     from lobio.klavinslab import sessions
#     registry = sessions.klregistry
#     registry.connector.api.set_timeout(120)
#     register_all(registry)
