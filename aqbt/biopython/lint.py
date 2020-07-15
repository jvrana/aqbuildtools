# TODO: This code is incomplete

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqFeature import SeqFeature, FeatureLocation
import inflection
import functools
import operator
from aqbt import bioadapter
from tqdm import tqdm
import hashlib
from typing import List, Tuple


def text_to_color(txt):
    hashed = hashlib.sha1(txt.encode('utf-8')).hexdigest()
    return "#" + hashed[-6:]


def parameterize(name):
    return inflection.parameterize(name)


def parts_df_to_records(df):
    SEQUENCE = "Sequence"
    SOURCE = "Source (Pubmed)"
    NAME = "Part Name"
    DESCRIPTION = "Part Description"
    ROLES = "SBOL_Roles"
    TYPES = "SBOL_Types"
    records = []
    for _, row in df.iterrows():
        seq_str = str(row[SEQUENCE]).strip()
        if seq_str and seq_str != 'nan':
            source = row[SOURCE]
            name = row[NAME]
            desc = row[DESCRIPTION]
            seq = Seq(row[SEQUENCE], generic_dna)
            roles = row[ROLES]
            types = row[TYPES]
            desc = "desc: {}, source: {}".format(desc, source)
            parametrized_name = parameterize(name)
            record = SeqRecord(seq, name=parametrized_name, id=parametrized_name, description=desc)
            info = {
                'source': source,
                'description': desc,
                'name_key': name
            }
            record.annotations.update(info)
            feature_loc = FeatureLocation(0, len(seq), strand=1)
            if isinstance(roles, list):
                role = roles[0]
            else:
                role = roles
            feature = SeqFeature(feature_loc, qualifiers=dict(info), type=role)
            feature.id = name
            feature.qualifiers.update({
                'ApEinfo_fwdcolor': [text_to_color(name)],
                'ApEinfo_revcolor': [text_to_color(name + "_r")],
                'ApEinfo_label': name,
                'label': name,
                'SBOLinfo_types': roles,
                'SBOLinfo_roles': types
            })
            feature.name = name
            record.features.append(feature)
            records.append(record)
        else:
            print("Skipping {}".format(row[NAME]))
    return records


def record_dict(records):
    by_key = {}
    for rec in records:
        key = rec.annotations['name_key'].lower()
        by_key.setdefault(key, list())
        if by_key[key]:
            raise ValueError("Key {} already exists".foramt(key))
        by_key[key] = rec
    return by_key


def record_from_parts(part_names, records):
    by_key = record_dict(records)
    try:
        new_record = functools.reduce(operator.add, [by_key[c.lower()] for c in part_names])
        return new_record
    except KeyError as e:
        raise e


def design_df_to_design_dict(melted_df):
    designs = {}

    for _, row in melted_df.iterrows():
        designs.setdefault(row['Design'], list())
        designs[row['Design']].append(row['Part'])
    return designs


def design_df_to_records(melted_df, parts):
    designs = design_df_to_design_dict(melted_df)

    new_records = []

    for design_name, part_names in designs.items():
        print(design_name)
        print(part_names)
        record = record_from_parts(part_names, parts)
        record.name = parameterize(design_name)
        new_records.append(record)

    return new_records


def update_benchling_existing(new_dnas, session, pbar=tqdm):
    def group_by_attr(dna_list, attr):
        by_dict = {}
        for dna in dna_list:
            key = getattr(dna, attr)
            by_dict.setdefault(key, list())
            by_dict[key].append(dna)
        return by_dict

    by_folder_id = group_by_attr(new_dnas, 'folder_id')

    for folder_id, new_dna_list in by_folder_id.items():
        existing_dna = session.DNASequence.list(folder_id=folder_id)
        existing_by_name = group_by_attr(existing_dna, 'name')

        iterator = new_dna_list
        if pbar:
            iterator = pbar(iterator)
            iterator.desc = "Adding to " + str(folder_id)

        for dna in iterator:
            if dna.name in existing_by_name:
                existing_dna = existing_by_name[dna.name][0]
                dna.id = existing_dna.id
                dna.update()
            else:
                dna.save()


import re
import bisect


class BioLinter(object):
    PROMOTER = ['promoter', 'http://identifiers.org/so/SO:0000167']
    TERMINATOR = ['terminator', 'http://identifiers.org/so/SO:0000141']

    def find_features(self, record: SeqRecord, feature_types: List[str]) -> List[SeqFeature]:
        features = []
        for f in record.features:
            if f.type.strip().lower() in [t.lower() for t in feature_types]:
                features.append(f)
        return features

    def find_promoter(self, record: SeqRecord) -> SeqFeature:
        return self.find_features(record, self.PROMOTER)

    def find_terminator(self, record: SeqRecord) -> SeqFeature:
        return self.find_features(record, self.TERMINATOR)

    def find_seq(self, record: SeqRecord, pattern: str, ignore_case: bool=False) -> List[re.Match]:
        if ignore_case:
            list(re.finditer(pattern, str(record.seq), re.IGNORECASE))
        return list(re.finditer(pattern, str(record.seq)))

    def find_start(self, record: SeqRecord) -> List[re.Match]:
        return self.find_seq(record, 'ATG', ignore_case=True)

    def find_stop(self, record: SeqRecord) -> List[re.Match]:
        return self.find_seq(record, 'TAG|TAA|TGA', ignore_case=True)

    @staticmethod
    def _by_frame(matches: List[re.Match]):
        frames = {
            0: [],
            1: [],
            2: []
        }
        for match in matches:
            frame = match.span()[0] % 3
            frames[frame].append((match.span()[0], match))

        for v in frames.values():
            v.sort()
        return frames

    def find_cds(self, record: SeqRecord, i: int, j: int) -> Tuple[int, int, int]:
        """
        Find the frame and location of the coding sequences. Returns a list of tuples of
        (frame {0, 1, 2}, start, end).

        :param record:
        :param i:
        :param j:
        :return:
        """
        region = record[i:j]
        starts = self.find_start(region)
        ends = self.find_stop(region)
        start_frames = self._by_frame(starts)
        end_frames = self._by_frame(ends)

        cds = []
        for f in [0, 1, 2]:
            _starts = start_frames[f]
            _ends = end_frames[f]

            for start in _starts:
                end_index = bisect.bisect_right([_x[0] for _x in _ends], start[0])
                if not end_index == len(_ends):
                    cds.append((f, start[0] + i, _ends[end_index][0] + i))
        return cds

    def validate_cds(self, record: SeqRecord, threshold: int=30) -> List[str]:
        """Validates a CDS. Ensures the start and stop codons are within the threshold of
        the promoter and terminator sequences."""
        errors = []
        promoters = self.find_promoter(record)
        terminators = self.find_terminator(record)
        if not promoters:
            errors.append('Error: {} has no promoters'.format(record.name))
        if not terminators:
            errors.append('Error: {} has no terminators'.format(record.name))

        if len(promoters) > 1:
            errors.append("Error: {} more than one promoter: {}".format(record.name, promoters))
        if len(terminators) > 1:
            errors.append("Error: {} more than one terminator: {}".format(record.name, terminators))
        if not errors:
            start_index = promoters[0].location.end
            end_index = terminators[0].location.start
            cds_list = self.find_cds(record, start_index, end_index)

            valid_cds = []
            for cds in cds_list:
                near_promoter = False
                near_terminator = False
                within_bp = threshold
                if start_index <= cds[1] and cds[1] <= start_index + within_bp:
                    near_promoter = True
                if end_index - within_bp <= cds[2] and cds[2] <= end_index:
                    near_terminator = True
                if (near_promoter and near_terminator):
                    valid_cds.append(cds)
            if not valid_cds:
                msg = "Error: No valid CDS"
                msg += "\n\tRecord: {}".format(record.name)
                msg += "\n\tPromoter End: {}".format(start_index)
                msg += "\n\tTerminator Start: {}".format(end_index)
                msg += "\n\tCDS List: {}".format((cds_list))
                errors.append(msg)

        return errors

    def errors(self, record):
        errors = []
        errors += self.validate_cds(record)
        return errors

    def lint(self, records):
        assert isinstance(records, list)
        for record in records:

            errs = self.errors(record)
            if errs:
                for err in errs:
                    print(err)
                print()


def prepare_for_gibson_assembly(record, parts):
    prefix = 'pp2'
    suffix = 'ts'
    promoter_suffix = 'ps'
    terminator_prefix = 'tp'

    parts_dict = record_dict(parts)

    linter = BioLinter()
    promoters = linter.find_promoter(record)
    terminators = linter.find_terminator(record)

    i = promoters[0].location.end
    j = terminators[0].location.start

    new_record = functools.reduce(operator.add, [
        parts_dict[prefix],
        record[:i],
        parts_dict[promoter_suffix],
        record[i:j],
        parts_dict[terminator_prefix],
        record[j:],
        parts_dict[suffix]
    ])
    return new_record