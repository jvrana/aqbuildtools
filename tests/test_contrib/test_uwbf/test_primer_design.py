from aqbt.contrib.uwbf import primers


def test_primer_df():
    template = primers.Templates.eyfp

    fwd_primers = [
        template[:20],
        template[:19],
        primers.rc(template[:30]),
        'agagataga',
        primers.rc(template[:30]),
    ]

    df = primers.create_anneal_df(template, primers=fwd_primers, names=list(range(5)))
    assert(len(df) == 4)


def test_primer_design():
    primer_design = primers.PrimerDesign()

    template = primers.Templates.eyfp

    pairs = primer_design.design_primers(template, fwd_primer='gtgagcaagggcgaggag')
    print(pairs)


def test_primer_design2():
    primer_design = primers.PrimerDesign()

    template = primers.Templates.eyfp

    pairs = primer_design.design_primers(template, start=200, length=200)
    print(pairs)
    assert pairs[0]['PAIR']['PRODUCT_SIZE'] == 200


def test_primer_design_from_list():
    primer_design = primers.PrimerDesign()

    template = primers.Templates.eyfp

    fwd_primers = [
        template[:20],
        template[:19],
        rc(template[:30])
    ]

    all_pairs = []
    for fwd in fwd_primers:
        pairs = primer_design.design_primers(template, fwd_primer=fwd)
        all_pairs += pairs
    print(all_pairs)


def test_anneal():
    template = primers.Templates.eyfp
    fwd_primers = [
        template[:20],
        template[:19],
        template[:30]
    ]
    df = primers.create_anneal_df(template, fwd_primers, list(range(3)))
    df = df[df['strand'] == 1]
    print(df)


def test_primer_anneal_mask():
    template = primers.Templates.eyfp
    fwd_primers = [
        template[:20],
        template[:19],
        template[:30],
        'agagaga'
    ]
    mask = primers.primer_anneal_mask(template, fwd_primers)
    print(mask)


def test_primer_anneal_mask():
    template = primers.Templates.eyfp
    fwd_primers = [
        template[:20],
        template[:19],
        primers.rc(template[:30]),
        'agagataga',
        primers.rc(template[:30]),
    ]

    mask, fwd_mask, rev_mask = primers.primer_anneal_mask(template, fwd_primers, ret_fwd_and_rev=True)
    assert mask == [True, True, True, False, True]
    assert fwd_mask == [True, True, False, False, False]
    assert rev_mask == [False, False, True, False, True]