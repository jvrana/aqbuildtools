from Bio.SeqRecord import SeqRecord
from pydent.models import Sample

from aqbt.aquarium.resolve_sequences import Resolver
import typing as t
import os
from Bio import SeqIO
import json
import pickle
import pytest
import dill

here = os.path.abspath(os.path.dirname(__file__))

def get_parents(yeast: Sample):
    parents = yeast.properties['Haploids']
    if not parents:
        parent = yeast.properties['Parent']
        if parent:
            parents.append(parent)
    return parents


def get_all_integrants(yeast: Sample, integrants: t.Optional[t.List[Sample]]=None):
    if integrants is None:
        integrants = list()
    integrant = yeast.properties.get('Integrant', None) or None
    if integrant:
        integrants.append(integrant)

    for p in get_parents(yeast):
        get_all_integrants(p, integrants)
    return integrants




def sample_iden(sample):
    return f"{sample.sample_type.name}:{sample.id}:{sample.name}"


def sanitize_genbank_name(x):
    import re

    x = re.sub('[<>\)\(]', '_', x)
    x = re.sub('\s', '_', x)
    x = x.strip('_')
    x = re.sub('[_]+', '_', x)
    return x

def test_sanitize_genbank_name():
    print(sanitize_genbank_name('URA-URA<pGRR-W10W8-URGR-W5> (v2)'))


@pytest.mark.parametrize('integrant_id',
                         [36483,
                          1380,
                          5535,
                          36332,
                          36353,
                          36482,
                          5533,
                          36471,
                          36339,
                          36475,
                          33580,
                          36451,
                          36355,
                          5534,
                          36445,
                          36450,
                          36446,
                          5536,
                          36447,
                          36449,
                          36468,
                          5537,
                          36474,
                          36470,
                          1096,
                          36469,
                          36473,
                          36472]
                         )
def test_resolve_sequences(integrant_id, aquarium, registry):
    filename = f'{integrant_id}.gb'
    if os.path.isfile(filename):
        print(f"{integrant_id} Done")
    else:
        resolver = Resolver(aquarium, registry)
        integrant = aquarium.Sample.find(integrant_id)
        seq = resolver.resolve_sequence(integrant)
        record: SeqRecord = resolver.registry.connector.convert(seq, frm="DNASequence", to="SeqRecord")
        if record is None:
            with open(filename + '.skipped', 'w') as f:
                ...
            raise RuntimeError()
        else:
            record.annotations['molecule_type'] = 'DNA'
            record.id = f'UWBF_{integrant_id}'
            record.name = f'UWBF_{integrant_id}'
            record.description = integrant.name
            with open(filename, 'w') as f:
                SeqIO.write([record], f, format='genbank')


@pytest.mark.parametrize('integrant_id',
                         [36425, 33889, 36264, 36241]
                         )
def test_resolve_sequence(integrant_id, aquarium, registry):
    filename = f'36474_{integrant_id}.gb'
    if os.path.isfile(filename):
        print(f"{integrant_id} Done")
    else:
        resolver = Resolver(aquarium, registry)
        integrant = aquarium.Sample.find(integrant_id)
        print(registry.find_in_registry(integrant))
        seq = resolver.resolve_sequence(integrant)
        record: SeqRecord = resolver.registry.connector.convert(seq, frm="DNASequence", to="SeqRecord")
        record.annotations['molecule_type'] = 'DNA'
        record.id = f'UWBF_{integrant_id}'
        record.name = f'UWBF_{integrant_id}'
        record.description = integrant.name
        with open(filename, 'w') as f:
            SeqIO.write([record], f, format='genbank')


@pytest.mark.parametrize('strain_id', [36564, 36565, 36566, 36567, 36568, 36569, 36570, 36571])
def test_resolve_all_integrants(strain_id, aquarium, registry):
    strain = aquarium.Sample.find(strain_id)
    integrants = get_all_integrants(strain)

    strain_json = {}
    strain_json['strains'] = []

    strain_json['strains'].append({
        'yeast': sample_iden(strain),
        'integrants': [sample_iden(x) for x in integrants]
    })
    with open(os.path.join(here, f'{strain_id}.strain.json'), 'w') as f:
        json.dump(strain_json, f, indent=2)

    integrants_skipped = []
    for integrant in integrants:
        integrant: Sample
        filepath = os.path.join(here, sample_iden(integrant) + '.pkl')
        if os.path.exists(filepath):
            print(f"file exists: {filepath}")
            continue
        resolver = Resolver(aquarium, registry)
        try:
            seq = resolver.resolve_sequence(integrant)
        except:
            integrants_skipped.append(integrant.id)
            continue
        if seq:
            record = resolver.registry.connector.convert(seq, frm="DNASequence", to="SeqRecord")
            _name = record.id or record.name
            record.name = sanitize_genbank_name(_name)
            record.id = sanitize_genbank_name(_name)
            with open(filepath, 'wb') as f:
                pickle.dump(record, f)
        else:
            integrants_skipped.append(integrant.id)
    with open(os.path.join(here, f'{sample_iden(strain)}.skipped.json'), 'w') as f:
        json.dump(integrants_skipped, f)