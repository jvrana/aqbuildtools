from pprint import pprint

from aqbt.contrib.uwbf import primer_utils


def test_primer_df():
    template = primer_utils.Templates.eyfp

    fwd_primers = [
        template[:20],
        template[:19],
        primer_utils.rc(template[:30]),
        "agagataga",
        primer_utils.rc(template[:30]),
    ]

    df = primer_utils.create_anneal_df(
        template, primers=fwd_primers, names=list(range(5))
    )
    assert len(df) == 4


def test_primer_design():
    primer_design = primer_utils.PrimerDesign()

    template = primer_utils.Templates.eyfp

    pairs = primer_design.design_cloning_primers(
        template, fwd_primer="gtgagcaagggcgaggag"
    )
    print(pairs)


def test_primer_design2():
    primer_design = primer_utils.PrimerDesign()

    template = primer_utils.Templates.eyfp

    pairs = primer_design.design_cloning_primers(template, start=200, length=200)
    print(pairs)
    assert pairs[0]["PAIR"]["PRODUCT_SIZE"] == 200


def test_anneal():
    template = primer_utils.Templates.eyfp
    fwd_primers = [template[:20], template[:19], template[:30]]
    df = primer_utils.create_anneal_df(template, fwd_primers, list(range(3)))
    df = df[df["strand"] == 1]
    print(df)


def test_primer_anneal_mask():
    template = primer_utils.Templates.eyfp
    fwd_primers = [template[:20], template[:19], template[:30], "agagaga"]
    mask = primer_utils.primer_anneal_mask(template, fwd_primers)
    print(mask)


def test_primer_anneal_mask():
    template = primer_utils.Templates.eyfp
    fwd_primers = [
        template[:20],
        template[:19],
        primer_utils.rc(template[:30]),
        "agagataga",
        primer_utils.rc(template[:30]),
    ]

    mask, fwd_mask, rev_mask = primer_utils.primer_anneal_mask(
        template, fwd_primers, ret_fwd_and_rev=True
    )
    assert mask == [True, True, True, False, True]
    assert fwd_mask == [True, True, False, False, False]
    assert rev_mask == [False, False, True, False, True]


def test_design_primers_from_list_1():
    template = primer_utils.Templates.eyfp
    primers = [
        "aaaaaaaaa",
        primer_utils.rc(template[-30:]),
    ]

    design = primer_utils.PrimerDesign()
    pairs = design.design_and_pick_primer(template, primers)
    assert not pairs[0]
    assert pairs[1]


def test_design_primers_from_list_2():
    template = primer_utils.Templates.eyfp
    primers = [
        "aaaaaaaaa",
        template[:30],
    ]

    design = primer_utils.PrimerDesign()
    pairs = design.design_and_pick_primer(template, primers)
    assert pairs[0]
    assert not pairs[1]


def test_design_primers_from_list_3():
    template = primer_utils.Templates.eyfp
    primers = ["aaaaaaaaa"]

    design = primer_utils.PrimerDesign()
    pairs = design.design_and_pick_primer(template, primers)
    assert not pairs[0]
    assert not pairs[1]


def test_design_fwd_primers_from_list():
    template = primer_utils.Templates.eyfp
    primers = [
        "aaaaaaaaa",
        primer_utils.rc(template[-30:]),
    ]

    design = primer_utils.PrimerDesign()
    pairs = design.design_fwd_and_pick_rev_primer(template, primers)
    assert pairs
    assert pairs[0]["RIGHT"]["META"]
    assert not pairs[0]["LEFT"]["META"]
    assert not design.design_rev_and_pick_fwd_primer(template, primers)


def test_design_rev_primers_from_list():
    template = primer_utils.Templates.eyfp
    primers = ["aaaaaaaaa", template[:30], "a"]

    design = primer_utils.PrimerDesign()
    pairs = design.design_rev_and_pick_fwd_primer(template, primers)
    assert pairs
    assert pairs[0]["LEFT"]["META"]
    assert not pairs[0]["RIGHT"]["META"]
    assert not design.design_fwd_and_pick_rev_primer(template, primers)


class TestAqPrimerDesign:
    def test_aq_primer_design__design_fwd_pick_rev(self, tools):
        aq = tools.sessions["default"]["aquarium"]
        df = primer_utils.get_aq_primers_df(aq, limit=100)

        design = primer_utils.PrimerDesign()
        template = primer_utils.Templates.eyfp + primer_utils.rc(list(df.sequence)[-1])

        pairs = design.design_fwd_and_pick_rev_primer(template, primer_seqs=df)
        assert pairs
        pprint(pairs)

    def test_aq_primer_design__design_rev_pick_fwd(self, tools):
        aq = tools.sessions["default"]["aquarium"]
        df = primer_utils.get_aq_primers_df(aq, limit=100)

        design = primer_utils.PrimerDesign()
        template = list(df.sequence)[-1] + primer_utils.Templates.eyfp

        pairs = design.design_rev_and_pick_fwd_primer(template, primer_seqs=df)
        assert pairs
        pprint(pairs)

    def test_aq_primer_design__design_fwd_pick_rev_with_overhang(self, tools):
        aq = tools.sessions["default"]["aquarium"]
        df = primer_utils.get_aq_primers_df(aq, limit=100)

        design = primer_utils.PrimerDesign()
        template = primer_utils.Templates.eyfp + primer_utils.rc(list(df.sequence)[-1])

        pairs = design.design_fwd_and_pick_rev_primer(
            template, primer_seqs=df, lflank="gaggattagata"
        )
        assert pairs
        assert pairs[0]["LEFT"]["OVERHANG"] == "gaggattagata"
        assert not pairs[0]["RIGHT"]["OVERHANG"]
        pprint(pairs)

    def test_aq_primer_design__design_rev_pick_fwd_with_overhang(self, tools):
        aq = tools.sessions["default"]["aquarium"]
        df = primer_utils.get_aq_primers_df(aq, limit=100)

        design = primer_utils.PrimerDesign()
        template = list(df.sequence)[-1] + primer_utils.Templates.eyfp
        pairs = design.design_rev_and_pick_fwd_primer(
            template, primer_seqs=df, rflank="gaggattagata"
        )
        assert pairs
        assert pairs[0]["RIGHT"]["OVERHANG"] == primer_utils.rc("gaggattagata")
        assert not pairs[0]["LEFT"]["OVERHANG"]
        pprint(pairs)
