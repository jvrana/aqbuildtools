from uuid import uuid4
import hashlib
from tqdm import tqdm_notebook as tqdm
import json

registry = sessions['default']['registry']
session = registry.session

PREFIX = "PlantTF_2020_Campaign_"
primer_type = session.SampleType.find_by_name('Primer')
fragment_type = session.SampleType.find_by_name('Fragment')
plasmid_type = session.SampleType.find_by_name('Plasmid')


def seq_sha1(seq: str) -> str:
    """Convert sequence string into a hash"""
    return hashlib.sha1(seq.strip().upper().encode()).hexdigest()


def new_name(molecule):
    name = molecule['__name__']
    if name == "PRIMER":
        typename = 'Primer'
    elif 'PCR' in name:
        typename = 'Fragment'
    elif name in ['GAP', 'SHARED_SYNTHESIZED_FRAGMENT']:
        typename = 'Synthesized'
    else:
        raise ValueError

    seqhash = seq_sha1(molecule['sequence']['bases'])
    return '{}_{}_{}'.format(typename, PREFIX, seqhash[-8:])


def new_sample(sample_type, name, description, project, properties):
    new_sample = sample_type.new_sample(
        name=name,
        description=description,
        project=project,
        properties=properties
    )
    #     sample = recycler.reuse(sample_type.id)
    sample = None
    if sample:
        for fv in sample.field_values:
            assert fv.value is None
            assert fv.child_sample_id is None
            assert fv.sample is None
        sample.name = new_sample.name
        sample.description = new_sample.description
        sample.update_properties(new_sample.properties)
        sample.needs_update = True
    else:
        sample = new_sample
    return sample


def _resolve_primer(m, inv_df):
    row = inv_df[inv_df['sequence_hash'] == seq_sha1(m['sequence']['bases'])]
    if len(row):
        sample_id = int(row['sample_id'].values[0])

        assert sample_id
        sample = session.Sample.find(sample_id)
        assert sample
        return sample
    else:
        return new_sample(sample_type=primer_type,
                          name=new_name(m),
                          description='dasi designed',
                          project='SD2',
                          properties={
                              'Anneal Sequence': m['__meta__']['SEQUENCE'],
                              'Overhang Sequence': m['__meta__']['OVERHANG'],
                              'T Anneal': round(m['__meta__']['TM'] - 2, 2)
                          }
                          )


def _resolve_fragment(m, results, df, r_dict, m_dict):
    reactions = m['used_as_output_to_reactions']
    assert len(reactions) == 1
    inputs = r_dict[reactions[0]]['inputs']
    input_molecules = [m_dict[i] for i in inputs]
    primers = [m for m in input_molecules if m['__name__'] == 'PRIMER']
    templates = [m for m in input_molecules if m['__name__'] == 'TEMPLATE']
    assert len(primers) == 2
    assert len(templates) == 1

    fwd = _resolve_primer(primers[0], df)
    rev = _resolve_primer(primers[1], df)

    return new_sample(
        sample_type=fragment_type,
        name=new_name(m),
        description='dasi designed',
        project='SD2',
        properties={
            'Forward Primer': fwd,
            'Reverse Primer': rev,
            'Template': _resolve_template(templates[0]),
            'Length': len(m['sequence']['bases']),
            'Sequence': '',
        }
    )


def _resolve_template(m):
    lims_id = m['sequence']['LIMS_ID']
    assert lims_id
    sample = session.Sample.find(lims_id)
    assert sample
    return sample


def _resolve_gblock(m):
    return new_sample(
        sample_type=fragment_type,
        name=new_name(m),
        description='dasi designed',
        project='SD2',
        properties={
            'Sequence': m['sequence']['bases'],
            'Length': len(m['sequence']['bases'])
        }
    )


def _resolve_plasmid(m):
    return new_sample(
        sample_type=plasmid_type,
        name=m['sequence']['name'],
        description='dasi designed',
        project='SD2',
        properties={
            'Sequence': m['sequence']['bases'],
            'Length': len(m['sequence']['bases']),
            'Bacterial Marker': 'Amp'
        }
    )


def _resolve_molecule(m, results, df, mol_dict, r_dict):
    if m['__name__'] in ['GAP', 'SHARED_SYNTHESIZED_FRAGMENT']:
        return _resolve_gblock(m)
    elif m['__name__'] in ['PCR_PRODUCT', 'PCR_PRODUCT_WITH_LEFT_PRIMER',
                           'PCR_PRODUCT_WITH_RIGHT_PRIMER', 'PCR_PRODUCT_WITH_PRIMERS']:
        return _resolve_fragment(m, results, df, r_dict, mol_dict)
    elif m['__name__'] == 'TEMPLATE':
        return _resolve_template(m)
    elif m['__name__'] == 'PRE-MADE DNA FRAGMENT':
        return _resolve_template(m)
    elif m['__name__'] == 'PRIMER':
        return _resolve_primer(m, df)
    elif m['__name__'] == 'PLASMID':
        return _resolve_plasmid(m)
    else:
        raise ValueError(m['__name__'] + " not recognized")


def resolve_molecules(results, df):
    molecules = results['molecules']
    mol_dict = {m['__index__']: m for m in results['molecules']}
    r_dict = {r['__index__']: r for r in results['reactions']}
    for m in molecules:
        m['__sample__'] = None

    for m in molecules:
        resolved = _resolve_molecule(m, results, df, mol_dict, r_dict)
        if resolved:
            m['__sample__'] = resolved
    return molecules