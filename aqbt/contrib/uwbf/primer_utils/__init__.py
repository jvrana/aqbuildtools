import functools
import re
from collections import namedtuple
from json import JSONDecodeError
from typing import Any
from typing import Dict
from typing import List
from typing import Optional
from typing import Tuple
from typing import Union

import pandas as pd
import primer3
import primer3plus
import pydent
from primer3plus.exceptions import Primer3PlusRunTimeError
from primer3plus.utils import reverse_complement as rc
from pydent import AqSession
from tqdm.auto import tqdm

bindings = namedtuple("Bindings", "fwd_picked rev_picked")
binding_mask = namedtuple("BindingMask", "mask fwd_mask rev_mask")


class PrimerException(Exception):
    ...


class Templates:
    eyfp = "gtgagcaagggcgaggagctgttcaccggggtggtgcccatcctggtcgagctggacggcgacgtaaacggccacaagttcagcgtgtccggcgagggcgagggcgatgccacctacggcaagctgaccctgaagttcatctgcaccaccggcaagctgcccgtgccctggcccaccctcgtgaccaccttcggctacggcctgcagtgcttcgcccgctaccccgaccacatgaagcagcacgacttcttcaagtccgccatgcccgaaggctacgtccaggagcgcaccatcttcttcaaggacgacggcaactacaagacccgcgccgaggtgaagttcgagggcgacaccctggtgaaccgcatcgagctgaagggcatcgacttcaaggaggacggcaacatcctggggcacaagctggagtacaactacaacagccacaacgtctatatcatggccgacaagcagaagaacggcatcaaggtgaacttcaagatccgccacaacatcgaggacggcagcgtgcagctcgccgaccactaccagcagaacacccccatcggcgacggccccgtgctgctgcccgacaaccactacctgagctaccagtccgccctgagcaaagaccccaacgagaagcgcgatcacatggtcctgctggagttcgtgaccgccgccgggatcactctcggcatggacgagctg"


def pagination(
    aqsession: AqSession,
    query,
    page_size,
    model="Sample",
    pbar=True,
    pbar_desc=None,
    limit=None,
):
    if limit:
        last_id = (
            aqsession.model_interface(model)
            .first(1, query=query, opts={"offset": limit})[0]
            .id
        )
    else:
        last_id = aqsession.model_interface(model).last(1, query=query)[0].id

    opts = {}
    if limit:
        opts["limit"] = limit
    iterator = aqsession.model_interface(model).pagination(
        query=query, page_size=page_size, opts=opts
    )
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
    for primers in pagination(
        sess,
        {"sample_type_id": sess.SampleType.find_by_name("Primer").id},
        page_size=page_size,
        limit=limit,
        pbar=pbar,
        pbar_desc="collecting primer seqs",
    ):
        sess.browser.get(
            primers, {"sample_type": "field_types", "field_values": "field_type"}
        )
        all_primers += primers
    return all_primers


def _clean_seq(seq):
    return re.sub(r"\s", "", seq)


def _primer_add_seqs(p):
    try:
        p.properties
        overhang = _clean_seq(p.properties["Overhang Sequence"] or "")
        anneal = _clean_seq(p.properties["Anneal Sequence"] or "")

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
        if hasattr(p, "anneal"):
            row = {
                "id": p.id,
                "name": p.name,
                "anneal": p.anneal,
                "overhang": p.overhang,
                "sequence": p.sequence,
            }
            primer_rows.append(row)
    return primer_rows


def create_primer_df(primers: List[pydent.models.Sample]) -> pd.DataFrame:
    """Create a pandas dataframe of the aquarium samples.

    :param primers:
    :return:
    """
    primer_rows = _primer_rows(primers)
    df = pd.DataFrame(primer_rows)
    df = df[df["sequence"] != ""]
    return df


def get_aq_primers_df(
    sess: AqSession, page_size=500, pbar=True, limit=None, timeout=60
) -> pd.DataFrame:
    primers = get_aq_primers(sess, page_size=page_size, pbar=pbar, limit=limit)
    return create_primer_df(primers)


def _is_left_end_terminal(df):
    return (df["strand"] == 1) & (df["start"] == 0)


def _is_right_end_terminal(df):
    template_len = len(df.meta["template"])
    return (df["strand"] == -1) & (df["start"] == template_len - 1)


def create_anneal_df(
    template: str,
    primers: List[str],
    names: List[str],
    n_bases: int = 12,
    min_tm=None,
    max_overhang=None,
    strand: Optional[int] = None,
    region: Optional[Tuple[int, int]] = None,
) -> pd.DataFrame:
    """Create a pandas Dataframe of the primer binding sites.

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
    assert len(primers) == len(names)
    meta = {
        "template": template,
        "n_bases": n_bases,
        "min_tm": min_tm,
        "max_overhang": max_overhang,
        "strand": strand,
        "region": region,
    }
    primer_list = list(zip(primers, names))
    fwd, rev = primer3plus.utils.anneal_iter(template, primer_list, n_bases=n_bases)

    rows = []
    for p in fwd + rev:
        seq = p["anneal"][-60:]
        tm = round(primer3.calcTm(seq), 2)
        p["tm"] = tm
        rows.append(p)

    df = pd.DataFrame(
        rows,
        columns=[
            "name",
            "anneal",
            "overhang",
            "primer",
            "start",
            "length",
            "top_strand_slice",
            "strand",
            "tm",
        ],
    )
    if strand is not None:
        assert strand in [1, -1]
        df = df[df["strand"] == strand]
    if min_tm is not None:
        df = df[df["tm"] >= min_tm]

    if region is not None:
        sel_a = df["top_strand_slice"].apply(lambda x: x[0] >= region[0])
        df = df[sel_a]
        sel_b = df["top_strand_slice"].apply(lambda x: x[1] < region[1])
        df = df[sel_b]
    if max_overhang is not None:
        sel = df["overhang"].apply(lambda x: len(x)) <= max_overhang
        df = df[sel]
    df.meta = meta
    df["left_term"] = _is_left_end_terminal(df)
    df["right_term"] = _is_right_end_terminal(df)
    return df


def _element_wise_and(a, b):
    return [_a and _b for _a, _b in zip(a, b)]


def primer_anneal_mask(
    template, primers, ret_fwd_and_rev: bool = False, **kwargs
) -> List[bool]:
    df = create_anneal_df(template, primers, list(range(len(primers))), **kwargs)
    mask = [False] * len(primers)
    fwd_mask = [False] * len(primers)
    rev_mask = fwd_mask[:]

    for _, row in df.iterrows():
        idx = row["name"]
        mask[idx] = True
        if row.strand == 1:
            fwd_mask[idx] = True
        elif row.strand == -1:
            rev_mask[idx] = True
    if ret_fwd_and_rev:
        return binding_mask(mask, fwd_mask, rev_mask)
    return mask


def _iter_mask_select(a, mask):
    assert len(a) == len(mask)
    for _a, _m in zip(a, mask):
        if mask:
            yield _a


def _mask_filter(a, mask):
    return list(_iter_mask_select(a, mask))


def parse_primer3_explain_flag(explain: dict):
    def parse_explain(explain):
        d = []
        for token in explain.split(", "):
            groups = re.search(r"(.+?)\s(\d+)", token).groups()
            d.append(groups)
        return dict(d)

    flag = {}
    for k, v in explain.items():
        if isinstance(v, str):
            v = parse_explain(v)
        flag[k] = v
    return flag


class PrimerDesign:
    default_params = dict(
        PRIMER_MAX_END_GC=3,
        PRIMER_MAX_POLY_X=4,
        PRIMER_WT_TM_LT=0.1,
        PRIMER_PAIR_WT_DIFF_TM=0.3,
        PRIMER_WT_HAIRPIN_TH=2.0,
        PRIMER_WT_SELF_ANY_TH=1.0,
        PRIMER_WT_SELF_END_TH=2.0,
        PRIMER_PAIR_WT_COMPL_ANY_TH=1.0,
        PRIMER_PAIR_WT_COMPL_END_TH=2.0,
    )

    def __init__(self, params: Optional[dict] = None):
        if params is None:
            params = self.default_params
        self.default_design_params = dict(params)

    @staticmethod
    def _resolve_region(template_len, start, length, end):
        """Resolve indices to a (start, len) tuple.

        :param template_len:
        :param start:
        :param length:
        :param end:
        :return:
        """
        if start is not None:
            if length is None and end is not None:
                length = end - start
                region = (start, length)
            elif length is not None:
                if end is not None:
                    assert end == start + length
                region = (start, length)
            else:
                raise ValueError
        else:
            region = (0, template_len)
        return region

    def new_cloning_primer_design(
        self,
        template,
        params=None,
        start=None,
        length=None,
        end=None,
        fwd_primer=None,
        rev_primer=None,
        lflank=None,
        rflank=None,
        min_anneal: int = 15,
    ):
        # if fwd_primer and lflank:
        #     raise ValueError
        # if rev_primer and rflank:
        #     raise ValueError
        if fwd_primer is not None:
            assert fwd_primer != ""
        if rev_primer is not None:
            assert rev_primer != ""
        start, product_length = self._resolve_region(len(template), start, length, end)
        region = (start, product_length)

        # primer design
        design = primer3plus.new(self.default_design_params)
        if params:
            design.update(params)
        design.settings.as_cloning_task()
        design.settings.template(template)
        design.settings.included((region[0], region[1]))
        design.settings.product_size((region[1], region[1]))
        if fwd_primer:
            design.settings.left_sequence(fwd_primer)
        if rev_primer:
            design.settings.right_sequence(rev_primer)
        if lflank:
            if fwd_primer:
                if not fwd_primer.startswith(lflank):
                    raise PrimerException("fwd_primer must start with lflank")
            left_overhang = lflank
        else:
            left_overhang = ""
        if rflank:
            if rev_primer:
                if not rev_primer.startswith(rc(rflank)):
                    raise PrimerException("rev_primer must start with rc(rflank)")
            right_overhang = rc(rflank)
        else:
            right_overhang = ""
        design.settings.left_overhang(left_overhang)
        design.settings.right_overhang(right_overhang)
        design.PRIMER_PICK_ANYWAY = False
        design.PRIMER_MIN_ANNEAL_CHECK = min_anneal
        design.settings.use_overhangs()
        design.settings.long_ok()
        return design

    def design_cloning_primers(
        self,
        template,
        start=None,
        length=None,
        end=None,
        fwd_primer=None,
        rev_primer=None,
        fwd_primer_meta=None,
        rev_primer_meta=None,
        lflank=None,
        rflank=None,
        min_anneal: int = 15,
        max_optimization_iterations: int = 3,
        n_return: int = 1,
    ) -> List[dict]:
        """Design primers for a specific region of a template. If no region
        information is provided, it will be assumed that the region includes
        the entire template.

        .. warning::

            Cyclic templates are not supported.

        :param template: template sequence
        :param start: optional start index of product (relative to template).
        :param length: optional length of product
        :param end: optional end index of product (relative to template)
        :param fwd_primer: (optional) Specific left primer sequence (5'->3') to use in design. The primer
            should bind to the template on the *bottom* strand; in other words, some portion of the right-hand side of the
            provided primer sequence should be found in the template sequence.
        :param rev_primer: (optional) Specific right primer sequence (5'->3') to use in design.
            The primer should bind to the template on the *top* strand; in other words, some portion of the
            right-hand side of the provided primer sequence should be found in the reverse complement of
            the template sequence. This is in the same orientation as you would provide to a synthesis vendor (e.g. IDT)
        :param fwd_primer_meta: additional meta data to attach to the forward primer ("LEFT")
        :param rev_primer_meta: additional meta data to attach to the reverse primer ("RIGHT")
        :param min_anneal: minimum number of bases to anneal
        :param max_optimization_iterations: maximum number of optimization iterations to perform
        :param n_return: number of primers to return
        :param lflank: (optional) For primer design, design primers such that this sequence flanks the left-hand
            side of the forward primer. If fwd_primer is also provided, the fwd_primer sequence *MUST* start with
            this lflank sequence.
        :param rflank: (optional) For primer design, design primers such that the product 3' end has this
            sequence. If rev_primer is also provided, the rev_primer sequence *MUST* start with the
            reverse complement of the rflank.
        :return:
        """
        design = self.new_cloning_primer_design(
            template=template,
            start=start,
            length=length,
            end=end,
            fwd_primer=fwd_primer,
            rev_primer=rev_primer,
            lflank=lflank,
            rflank=rflank,
            min_anneal=min_anneal,
        )
        design.settings.primer_num_return(n_return)

        pairs, explain = design.run_and_optimize(
            max_iterations=max_optimization_iterations, pick_anyway=True
        )
        pairs = list(pairs.values())
        results = []
        fwd_primer_meta = fwd_primer_meta or {}
        rev_primer_meta = rev_primer_meta or {}
        for x in pairs:
            x["LEFT"]["META"] = dict(fwd_primer_meta)
            x["RIGHT"]["META"] = dict(rev_primer_meta)
            results.append(x)
        # explain = parse_primer3_explain_flag(explain)
        return results

    def design_and_pick_primer(
        self,
        template,
        primer_seqs,
        primer_meta_list=None,
        start=None,
        end=None,
        length=None,
        **kwargs
    ):
        if isinstance(primer_seqs, pd.DataFrame):
            df = primer_seqs
            primer_seqs = df.sequence
            primer_meta_list = df.T.to_dict().values()
        if primer_meta_list is None:
            primer_meta_list = [{"index": i} for i in list(range(len(primer_seqs)))]
        r = self._resolve_region(len(template), start=start, end=end, length=length)
        _, fwd_mask, rev_mask = primer_anneal_mask(
            template[r[0] : r[0] + r[1]], primer_seqs, ret_fwd_and_rev=True
        )

        all_pairs = bindings([], [])
        for i, (primer_key, meta_key, mask) in enumerate(
            [
                ("fwd_primer", "fwd_primer_meta", fwd_mask),
                ("rev_primer", "rev_primer_meta", rev_mask),
            ]
        ):
            mfilter = functools.partial(_mask_filter, mask=mask)
            _kwargs = dict(kwargs)
            for seq, meta in zip(mfilter(primer_seqs), mfilter(primer_meta_list)):
                _kwargs[primer_key] = seq
                _kwargs[meta_key] = meta
                try:
                    pairs = self.design_cloning_primers(
                        template, start=start, end=end, length=length, **_kwargs
                    )
                    all_pairs[i].extend(pairs)
                except Primer3PlusRunTimeError:
                    pass
                except PrimerException:
                    pass
        all_pairs[0].sort(key=lambda x: x["PAIR"]["PENALTY"])
        all_pairs[1].sort(key=lambda x: x["PAIR"]["PENALTY"])
        return all_pairs

    def design_fwd_and_pick_rev_primer(
        self,
        template,
        primer_seqs,
        primer_meta_list=None,
        start=None,
        end=None,
        length=None,
        **kwargs
    ):
        result = self.design_and_pick_primer(
            template,
            primer_seqs,
            primer_meta_list,
            start=start,
            end=end,
            length=length,
            **kwargs
        )
        return result.rev_picked

    def design_rev_and_pick_fwd_primer(
        self,
        template,
        primer_seqs,
        primer_meta_list=None,
        start=None,
        end=None,
        length=None,
        **kwargs
    ):
        result = self.design_and_pick_primer(
            template,
            primer_seqs,
            primer_meta_list,
            start=start,
            end=end,
            length=length,
            **kwargs
        )
        return result.fwd_picked


class AqPrimerDesign(PrimerDesign):
    def __init__(self, session, params=None):
        super().__init__(params)
        self.session = session
        self._primer_df = None

    @property
    def primer_df(self):
        if self._primer_df is None:
            self._primer_df = get_aq_primers_df(self.session)
        return self._primer_df

    def reset(self):
        self._primer_df = None

    def pick_fwd_aq_primer(self):
        self.design_cloning_primers(
            template, self.primer_df.sequence, self.primer_df.name
        )
