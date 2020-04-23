import itertools
from typing import List

from primer3plus.design import primer3
from pydent.models import Sample

from aqbt import biopython


class LintError:

    WARNING = "warning"
    ERROR = "error"
    TYPES = [WARNING, ERROR]

    def __init__(self, err_type: str, sample: Sample, msg: str):
        self.msg = msg
        self.sample = sample
        self.err_type = err_type
        assert self.err_type in self.TYPES

    @classmethod
    def warning(cls, sample: Sample, msg: str):
        return cls(cls.WARNING, sample, msg)

    @classmethod
    def error(cls, sample: Sample, msg: str):
        return cls(cls.ERROR, sample, msg)

    def __str__(self):
        return "{err} {sample_id}:{sample_name}: {msg}".format(
            err=self.err_type.upper(),
            sample_id=self.sample.id,
            sample_name=self.sample.name,
            msg=self.msg,
        )


class Linter:
    def __init__(self):
        self.PRIMER_TM_FWD_REV_DIFF_THRESHOLD = 5
        self.PRIMER_TA_TO_TM = 3
        self.PRIMER_TA_DIFF_THRESHOLD = 3
        self.LOW_TA_THRESHOLD = 50
        self.LENGTH_DIFF_RATIO_THRESHOLD = 0.05
        self.errors = []

    def warning(self, sample: Sample, msg: str):
        self.errors.append(LintError.warning(sample, msg))

    def error(self, sample: Sample, msg: str):
        self.errors.append(LintError.error(sample, msg))

    @staticmethod
    def _wrap_messages(msg_type: str, prefix: str, sample, msgs: List[str]):
        return [
            "{msg_type} {prefix}{sample}: {msg}".format(
                msg_type=msg_type,
                prefix=prefix,
                sample="{}:{}".format(sample.id, sample.name[:20]),
                msg=msg,
            )
            for msg in msgs
        ]

    def lint_primer(self, sample):
        anneal = sample.properties["Anneal Sequence"] or ""
        # overhang = sample.properties["Overhang Sequence"] or ""
        ta = sample.properties["T Anneal"]

        errors = []

        if not anneal:
            errors.append("Annealing sequence is missing")

        tm_anneal = primer3.calcTm(anneal)

        calc_ta = tm_anneal - self.PRIMER_TA_TO_TM
        if abs(calc_ta - ta) > self.PRIMER_TA_DIFF_THRESHOLD:
            self.warning(
                sample,
                "Reported 'T Anneal' {} does not match calculated T Anneal {}".format(
                    ta, calc_ta
                ),
            )

        if calc_ta < self.LOW_TA_THRESHOLD:
            self.warning(
                sample,
                "Calculated 'T Anneal' {} lower that threshold {}".format(
                    calc_ta, self.LOW_TA_THRESHOLD
                ),
            )

        return errors

    def lint_fragment(self, registry, fragment):
        errors = []

        fwd = registry.get_primer_sequence(fragment.properties["Forward Primer"])
        rev = registry.get_primer_sequence(fragment.properties["Reverse Primer"])
        template = registry.get_sequence(fragment.properties["Template"])
        reported_sequence = registry.get_sequence(fragment)
        if not reported_sequence:
            self.warning(fragment, "There is no reported sequence")

        if not fwd:
            self.error(
                fragment,
                "Forward Primer '{}' missing sequence".format(
                    fragment.properties["Forward Primer"].name
                ),
            )
        if not rev:
            self.error(
                fragment,
                "Reverse Primer '{}' missing sequence".format(
                    fragment.properties["Reverse Primer"].name
                ),
            )
        if not template:
            self.error(
                fragment,
                "Template '{}' is missing sequence".format(
                    fragment.properties["Template"].name
                ),
            )

        if not template or not fwd or not rev:
            return []

        self.lint_primer(fragment.properties["Forward Primer"])
        self.lint_primer(fragment.properties["Reverse Primer"])

        template_record = registry.connector.convert(template, to="SeqRecord")

        products, fwd_matches, rev_matches = biopython.pcr_amplify(
            (fwd, rev),
            template_record,
            cyclic=template.is_circular,
            name=fragment.name,
            return_matches=True,
        )

        if not fwd_matches:
            self.error(
                fragment,
                "There are no forward matches with '{}' and '{}'".format(
                    fragment.properties["Forward Primer"].name,
                    fragment.properties["Reverse Primer"].name,
                ),
            )
        if not rev_matches:
            self.error(
                fragment,
                "There are no reverse anneal with '{}' and '{}'".format(
                    fragment.properties["Forward Primer"].name,
                    fragment.properties["Reverse Primer"].name,
                ),
            )
        reported_fwd_ta = fragment.properties["Forward Primer"].properties["T Anneal"]
        reported_rev_ta = fragment.properties["Reverse Primer"].properties["T Anneal"]

        for f, r in itertools.product(fwd_matches, rev_matches):
            fwd_anneal = f["anneal"]
            rev_anneal = r["anneal"]
            fwd_ta = primer3.calcTm(fwd_anneal[-60:]) - self.PRIMER_TA_TO_TM
            rev_ta = primer3.calcTm(rev_anneal[-60:]) - self.PRIMER_TA_TO_TM
            if abs(fwd_ta - rev_ta) > self.PRIMER_TM_FWD_REV_DIFF_THRESHOLD:
                errors.append(
                    LintError.warning(
                        fragment,
                        "Forward and reverse primer have different "
                        "annealing temperatures: {} vs {}".format(fwd_ta, rev_ta),
                    )
                )
            if abs(reported_fwd_ta - fwd_ta) > self.PRIMER_TA_DIFF_THRESHOLD:
                errors.append(
                    LintError.warning(
                        fragment.properties["Forward Primer"],
                        "Reported 'T Anneal' {} does not match calculated T Anneal {}"
                        " for the specified template".format(reported_fwd_ta, fwd_ta),
                    )
                )
            if abs(reported_rev_ta - rev_ta) > self.PRIMER_TA_DIFF_THRESHOLD:
                errors.append(
                    LintError.warning(
                        fragment.properties["Reverse Primer"],
                        "Reported 'T Anneal' {} does not match calculated T Anneal {}"
                        " for the specified template".format(reported_rev_ta, rev_ta),
                    )
                )

        if len(products) > 1:
            self.warning(
                fragment,
                "There is more than one product. Found {} products.".format(
                    len(products)
                ),
            )
        elif len(products) == 0:
            self.error(fragment, "There are no products")
            return

        intended_product = products[0][0]

        if reported_sequence:
            if (
                str(intended_product.seq).upper()
                != str(reported_sequence.bases).upper()
            ):
                self.error(
                    fragment,
                    "There is a conflict between the expected and reported"
                    " sequence.",
                )

        reported_length = fragment.properties["Length"]
        expected_length = len(intended_product.seq)
        diff_length = abs(reported_length - expected_length)
        if diff_length / reported_length > self.LENGTH_DIFF_RATIO_THRESHOLD:
            self.error(
                fragment,
                "The reported length {} is more than"
                " {}% different the the expected length {}".format(
                    reported_length, self.LENGTH_DIFF_RATIO_THRESHOLD, expected_length
                ),
            )
        elif diff_length:
            self.warning(
                fragment,
                "The reported length {} is different than the"
                " expected length {}".format(reported_length, expected_length),
            )

        return products

    def report(self):
        warnings = [e for e in self.errors if e.err_type == LintError.WARNING]
        errors = [e for e in self.errors if e.err_type == LintError.ERROR]
        warnings = list(set(warnings))
        errors = list(set(errors))
        for e in errors:
            print(e)
        for w in warnings:
            print(w)
