# ----------------------------------------------------------------------------
# Copyright (c) 2026, Oddant1.
#
# No license has been chosen for this software yet.
#
# See the LICENSE file distributed with this software.
# ----------------------------------------------------------------------------

"""Allele calling adapted from Real_AlleleV11 by Jason Sahl and Grant
Pemberton.

Sequence IDs are expected to look like ``<sample>_<read count>``: everything
before the last underscore identifies the sample, and the final field is the
number of reads supporting that consensus sequence. Alternatively, a feature
table can give the read count of each sequence in each sample, in which case
the sequence IDs are only feature IDs.
"""

import collections
import re
from typing import NamedTuple

import pandas as pd
import rachis


_READ_COUNT = re.compile(r'[0-9]+')
_ALLELE_ID = re.compile(r'A([0-9]+)')


class _Record(NamedTuple):
    id: str
    sample: str
    read_count: int
    sequence: str


def call_alleles(sequences: pd.Series,
                 reference_alleles: pd.Series = None,
                 table: pd.DataFrame = None,
                 min_coverage: int = 10,
                 proportion: float = 0.5) -> (pd.Series,
                                              pd.Series,
                                              rachis.Metadata,
                                              rachis.Metadata,
                                              rachis.Metadata):
    if table is None:
        records = _parse_records(sequences)
    else:
        records = _records_from_table(table, sequences)
    passing = _filter_by_proportion(
        _filter_by_coverage(records, min_coverage), proportion)
    if not passing:
        raise ValueError(
            'None of the %d consensus sequences has at least %d reads '
            '(min_coverage), so there is nothing to genotype.'
            % (len(records), min_coverage))

    frequencies = _count_alleles(passing)
    names = _assign_allele_names(
        frequencies, _load_reference_alleles(reference_alleles))

    # Real_AlleleV11 reports the passing records in input order. Records are
    # matched whole rather than by ID, because with a table two of a
    # sample's sequences can have the same read count, and so the same ID.
    kept = set(passing)
    passing = [record for record in records if record in kept]

    alleles = pd.Series({allele: sequence
                         for sequence, allele in names.items()})
    genotyped_sequences = pd.Series(
        {'%s_%s' % (record.id, names[record.sequence]): record.sequence
         for record in passing})
    genotypes = _tabulate_genotypes(passing, names)
    allele_frequencies = pd.DataFrame(
        {'frequency': [count for _, count in frequencies]},
        index=pd.Index([names[sequence] for sequence, _ in frequencies],
                       name='feature-id'))

    return (alleles,
            genotyped_sequences,
            rachis.Metadata(genotypes),
            rachis.Metadata(allele_frequencies),
            rachis.Metadata(_summarize_zygosity(genotypes)))


def _parse_records(sequences):
    records = []
    malformed = []
    for id_, sequence in sequences.items():
        # Real_AlleleV11 never keeps consensus records, whatever their count.
        if 'CONSENSUS' in id_:
            continue
        sample, _, read_count = id_.rpartition('_')
        if not sample or not _READ_COUNT.fullmatch(read_count):
            malformed.append(id_)
            continue
        records.append(
            _Record(id_, sample, int(read_count), str(sequence)))

    if malformed:
        raise ValueError(
            'Sequence IDs must have the form <sample>_<read count>, where '
            'the read count is a whole number (e.g. SampleA_S1_1500). These '
            'IDs do not: %s' % _format_ids(malformed))
    return records


def _records_from_table(table, sequences):
    missing = [feature for feature in table.columns
               if feature not in sequences.index]
    if missing:
        raise ValueError(
            'Every feature in the table must be one of the sequences. These '
            'features are not: %s' % _format_ids(missing))

    records = []
    for sample, counts in table.iterrows():
        for feature, count in counts[counts > 0].items():
            read_count = int(count)
            # Named like a consensus sequence ID, so that the genotyped
            # sequences are named the same way whichever input was used.
            records.append(_Record('%s_%d' % (sample, read_count), sample,
                                   read_count, str(sequences[feature])))
    return records


def _filter_by_coverage(records, min_coverage):
    return [record for record in records if record.read_count >= min_coverage]


def _filter_by_proportion(records, proportion):
    by_sample = {}
    for record in records:
        by_sample.setdefault(record.sample, []).append(record)

    # Keep each sample's most abundant sequence, plus the runner-up when it
    # has at least `proportion` as many reads. Any others are discarded.
    passing = []
    for sample_records in by_sample.values():
        top, *others = sorted(sample_records,
                              key=lambda record: record.read_count,
                              reverse=True)
        passing.append(top)
        if others and others[0].read_count / top.read_count >= proportion:
            passing.append(others[0])
    return passing


def _count_alleles(records):
    # Most frequent first; ties stay in the order they were first seen.
    counts = collections.Counter(record.sequence for record in records)
    return sorted(counts.items(), key=lambda item: item[1], reverse=True)


def _load_reference_alleles(reference_alleles):
    if reference_alleles is None:
        return {}

    known = {}
    bad_ids = []
    duplicates = []
    for allele, sequence in reference_alleles.items():
        sequence = str(sequence)
        if not _ALLELE_ID.fullmatch(allele):
            bad_ids.append(allele)
        elif sequence in known:
            duplicates.append('%s and %s' % (known[sequence], allele))
        else:
            known[sequence] = allele

    if bad_ids:
        raise ValueError(
            'Reference allele IDs must have the form A<number> (e.g. A12), '
            'so that new alleles can be numbered after them. These IDs do '
            'not: %s' % _format_ids(bad_ids))
    if duplicates:
        raise ValueError(
            'Each reference allele must have a unique sequence. These '
            'reference alleles share a sequence: %s'
            % _format_ids(duplicates, quote=False))
    return known


def _assign_allele_names(frequencies, known):
    """Map each sequence to an allele name.

    Sequences in `known` keep their reference names; new sequences are named
    in order of decreasing frequency, starting after the highest reference
    number. The result lists the reference alleles first, so it can serve as
    the reference for a later run.
    """
    names = dict(known)
    next_number = 1 + max(
        (int(_ALLELE_ID.fullmatch(allele).group(1))
         for allele in known.values()), default=0)
    for sequence, _ in frequencies:
        if sequence not in names:
            names[sequence] = 'A%d' % next_number
            next_number += 1
    return names


def _tabulate_genotypes(records, names):
    by_sample = {}
    for record in records:
        by_sample.setdefault(record.sample, []).append(record)

    rows = []
    for sample, sample_records in by_sample.items():
        first, *second = sorted(sample_records,
                                key=lambda record: record.read_count,
                                reverse=True)
        # As in Real_AlleleV11, homozygotes list their only allele twice.
        second = second[0] if second else first
        rows.append((sample,
                     first.read_count, names[first.sequence],
                     second.read_count, names[second.sequence]))

    return pd.DataFrame(
        rows, columns=['sample-id', 'read-count-1', 'allele-1',
                       'read-count-2', 'allele-2']).set_index('sample-id')


def _summarize_zygosity(genotypes):
    total = len(genotypes)
    heterozygous = int((genotypes['allele-1'] != genotypes['allele-2']).sum())
    counts = [total - heterozygous, heterozygous]
    return pd.DataFrame(
        {'sample-count': counts,
         'percentage': [round(100 * count / total, 2) for count in counts]},
        index=pd.Index(['homozygous', 'heterozygous'], name='id'))


def _format_ids(ids, limit=5, quote=True):
    shown = [repr(id_) if quote else id_ for id_ in ids[:limit]]
    if len(ids) > limit:
        shown.append('and %d more' % (len(ids) - limit))
    return ', '.join(shown)
