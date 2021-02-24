from typing import List, Optional, Tuple, Union, Dict, Any

from tqdm.auto import tqdm
import re
from json import JSONDecodeError
import pandas as pd
import primer3
import primer3plus
import pydent
from pydent import AqSession


def pagination(aqsession: AqSession, query, page_size, model="Sample", pbar=True, pbar_desc=None, limit=None):
    if limit:
        last_id = aqsession.model_interface(model).first(1, query=query, opts={'offset': limit})[0].id
    else:
        last_id = aqsession.model_interface(model).last(1, query=query)[0].id

    opts = {}
    if limit:
        opts['limit'] = limit
    iterator = aqsession.model_interface(model).pagination(query=query, page_size=page_size, opts=opts)
    if pbar:
        pbar_iter = tqdm(iterator, total=last_id)

        if pbar_desc:
            pbar_iter.set_description(pbar_desc)
        last_max_id = -1
        for primers in pbar_iter:
            ids = [m.id for m in primers]
            max_id = max(ids)
            pbar_iter.update(max_id - last_max_id)
            last_max_id = max_id
            yield primers
    else:
        yield from iterator


def get_aq_primers(sess, page_size=500, pbar=True, limit=None):
    all_primers = []
    for primers in pagination(sess, {'sample_type_id': sess.SampleType.find_by_name('Primer').id}, page_size=page_size,
                              limit=limit, pbar=pbar, pbar_desc='collecting primer seqs'):
        sess.browser.get(primers, {'sample_type': 'field_types', 'field_values': 'field_type'})
        all_primers += primers
    return all_primers


def _clean_seq(seq):
    return re.sub('\s', '', seq)


def _primer_add_seqs(p):
    try:
        p.properties
        overhang = _clean_seq(p.properties['Overhang Sequence'] or '')
        anneal = _clean_seq(p.properties['Anneal Sequence'] or '')

        sequence = overhang + anneal
        p.overhang = overhang
        p.anneal = anneal
        p.sequence = sequence
    except JSONDecodeError as e:
        pass


def _primer_rows(primers):
    for p in primers:
        _primer_add_seqs(p)

    primer_rows = []
    for p in primers:
        if hasattr(p, 'anneal'):
            row = {
                'id': p.id,
                'name': p.name,
                'anneal': p.anneal,
                'overhang': p.overhang,
                'sequence': p.sequence
            }
            primer_rows.append(row)
    return primer_rows


def create_primer_df(primers: List[pydent.models.Sample]) -> pd.DataFrame:
    """
    Create a pandas dataframe of the aquarium samples.

    :param primers:
    :return:
    """
    primer_rows = _primer_rows(primers)
    df = pd.DataFrame(primer_rows)
    return df


def get_aq_primers_df(sess: AqSession,  page_size=500, pbar=True, limit=None, timeout=60) -> pd.DataFrame:
    primers = get_aq_primers(sess, page_size=page_size, pbar=pbar, limit=limit)
    return create_primer_df(primers)


def _is_left_end_terminal(df):
    return (df['strand'] == 1) & (df['start'] == 0)


def _is_right_end_terminal(df):
    template_len = len(df.meta['template'])
    return (df['strand'] == -1) & (df['start'] == template_len - 1)


def create_anneal_df(template: str, primers: List[str], names: List[str], n_bases: int=12, min_tm=None,
              max_overhang=None, strand: Optional[int]=None, region: Optional[Tuple[int, int]]=None) -> pd.DataFrame:
    """
    Create a pandas Dataframe of the primer binding sites.

    :param template: tempalte sequence
    :param primers: primer sequence list
    :param names: primer name list
    :param n_bases: minimum number of bases to anneal
    :param min_tm: minimum tm
    :param max_overhang: maximum overhang length
    :param strand: downselect the strand (1, -1)
    :param region: downselect the region [a, b)
    :return:
    """
    meta = {
        'template': template,
        'n_bases': n_bases,
        'min_tm': min_tm,
        'max_overhang': max_overhang,
        'strand': strand,
        'region': region
    }
    primer_list = list(zip(primers, names))
    fwd, rev = primer3plus.utils.anneal_iter(template, primer_list, n_bases=n_bases)

    filtered_fwd, filtered_rev = [], []

    rows = []
    for p in fwd + rev:
        seq = p['anneal'][-60:]
        tm = round(primer3.calcTm(seq), 2)
        p['tm'] = tm
        rows.append(p)

    df = pd.DataFrame(rows)
    if strand is not None:
        assert strand in [1, -1]
        df = df[df['strand'] == strand]
    if min_tm is not None:
        df = df[df['tm'] >= min_tm]

    if region is not None:
        sel_a = df['top_strand_slice'].apply(lambda x: x[0] >= region[0])
        df = df[sel_a]
        sel_b = df['top_strand_slice'].apply(lambda x: x[1] < region[1])
        df = df[sel_b]
    if max_overhang is not None:
        sel = df['overhang'].apply(lambda x: len(x)) <= max_overhang
        df = df[sel]
    df.meta = meta
    df['left_term'] = _is_left_end_terminal(df)
    df['right_term'] = _is_right_end_terminal(df)
    return df


def _design_primers(
        template: str,
        start: int,
        length: int,
        lseq: Union[None, str],
        rseq: Union[None, str],
        left_overhang: Union[None, str] = None,
        right_overhang: Union[None, str] = None,
) -> Tuple[Dict[int, dict], Dict[str, Any]]:
    """Design primers flanking the specified.

    :class:`Region.<dasi.utils.Region>`. If the region is cyclic and spans the
    origin, this method will handle the appropriate manipulations to design
    primers around the origin and restore the locations of the resulting primer
    pairs.

    :param template: the template string to design primers
    :param region: region specified to design primers around. Regions are exclusive at
                    their end points (`.b` parameter)
    :param lseq: optionally provided left sequence
    :param rseq: optionally provided right sequence
    :param left_overhang: optionally provided left overhang sequence of the primer
    :param right_overhang: optionally provided right overhang sequence of the primer
    :return: tuple of pairs and the 'explain' dictionary.
    """
    design = primer3plus.new()
    design.settings.as_cloning_task()

    #     if region.direction == -1:
    #         region = region.flip()
    #         template = rc(template)

    if lseq and left_overhang:
        raise ValueError
    if rseq and right_overhang:
        raise ValueError

    region = (start, length)

    spans_origin = False
    if spans_origin:
        ...
        # adjusted_template = region.get_slice(template) + region.invert()[0].get_slice(
        #     template
        # )
        # design.settings.template(adjusted_template)
        # design.settings.included((0,))
        # index = list(region) + list(region.invert()[0])
    else:
        design.settings.template(template)
        design.settings.included((region[0], region[1]))
        index = None

    if lseq:
        design.settings.left_sequence(lseq)
    if rseq:
        design.settings.right_sequence(rseq)

    if left_overhang is None:
        left_overhang = ""
    if right_overhang is None:
        right_overhang = ""

    design.settings.product_size((region[1], region[1]))
    design.settings.left_overhang(left_overhang)
    design.settings.right_overhang(right_overhang)
    design.PRIMER_PICK_ANYWAY = False
    design.PRIMER_MIN_ANNEAL_CHECK = 15
    design.settings.use_overhangs()
    design.settings.long_ok()

    design.logger.set_level("INFO")

    # TODO: remove debugging code
    try:
        pairs, explain = design.run_and_optimize(
            max_iterations=3, pick_anyway=True
        )
    except Exception as e:
        import json

        raise e
    if index is not None:
        for pair in pairs.values():
            loc = pair["LEFT"]["location"]
            pair["LEFT"]["location"] = (index[loc[0]], loc[1])

            loc = pair["RIGHT"]["location"]
            pair["RIGHT"]["location"] = (index[loc[0]], loc[1])
    return pairs, explain