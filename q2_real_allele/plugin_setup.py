# ----------------------------------------------------------------------------
# Copyright (c) 2026, Oddant1.
#
# No license has been chosen for this software yet.
#
# See the LICENSE file distributed with this software.
# ----------------------------------------------------------------------------

from rachis.plugin import Citations, Float, Int, Plugin, Range
from q2_types.feature_data import FeatureData, Sequence
from q2_types.feature_table import FeatureTable, Frequency
from q2_types.metadata import ImmutableMetadata

from q2_real_allele import __version__
from q2_real_allele._methods import call_alleles
from q2_real_allele._examples import (call_alleles_basic,
                                      call_alleles_with_reference)

citations = Citations.load("citations.bib", package="q2_real_allele")

plugin = Plugin(
    name="real-allele",
    version=__version__,
    website="https://example.com",
    package="q2_real_allele",
    description="Call alleles and genotypes from per-sample consensus "
                "amplicon sequences using the Real Allele algorithm.",
    short_description="Allele and genotype calling from consensus sequences.",
    citations=[citations['Stone2024']]
)

plugin.methods.register_function(
    function=call_alleles,
    inputs={'sequences': FeatureData[Sequence],
            'reference_alleles': FeatureData[Sequence],
            'table': FeatureTable[Frequency]},
    parameters={'min_coverage': Int % Range(1, None),
                'proportion': Float % Range(0, 1, inclusive_end=True)},
    outputs=[('alleles', FeatureData[Sequence]),
             ('genotyped_sequences', FeatureData[Sequence]),
             ('genotype_table', ImmutableMetadata),
             ('allele_frequencies', ImmutableMetadata),
             ('heterozygosity', ImmutableMetadata)],
    input_descriptions={
        'sequences': (
            'Consensus sequences to genotype. Unless a `table` is provided, '
            'each ID must have the form <sample>_<read count>: everything '
            'before the last underscore identifies the sample, and the final '
            'field is the number of reads supporting the sequence (e.g. '
            'RHIMI-ARA13-xx-DN-020-JB_S86_17441), and sequences whose ID '
            'contains "CONSENSUS" are ignored. With a `table`, these are the '
            'features of the table, and their IDs can be anything.'),
        'reference_alleles': (
            'Previously named alleles, such as the `alleles` output of an '
            'earlier run. Sequences identical to a reference allele are '
            'given its name, and new alleles are numbered after the highest '
            'reference allele. IDs must have the form A<number> (e.g. A1).'),
        'table': (
            'The number of reads supporting each of the `sequences` in each '
            'sample, with the sequences as features. When provided, samples '
            'and read counts are taken from this table instead of from the '
            'sequence IDs, which also allows two sequences in a sample to '
            'have the same read count. Sequences with a count of zero in a '
            'sample are not part of that sample.')
    },
    parameter_descriptions={
        'min_coverage': (
            'Consensus sequences supported by fewer reads than this are '
            'discarded.'),
        'proportion': (
            'Only the two most abundant consensus sequences of each sample '
            'are considered. The second is kept, making the sample '
            'heterozygous, if its read count is at least this proportion of '
            'the first\'s; otherwise the sample is homozygous for the most '
            'abundant sequence.')
    },
    output_descriptions={
        'alleles': (
            'Allele sequences, named A1, A2, ... in order of decreasing '
            'frequency. If reference alleles were provided, this contains '
            'all of them followed by any new alleles, so it can be used as '
            'the reference for the next run.'),
        'genotyped_sequences': (
            'The consensus sequences that passed filtering, with the name of '
            'their allele appended to the ID (e.g. <sample>_<read '
            'count>_A1).'),
        'genotype_table': (
            'The genotype of each sample: the read count and allele of its '
            'most abundant (1) and second (2) consensus sequence. '
            'Homozygous samples list the same read count and allele twice.'),
        'allele_frequencies': (
            'The number of retained consensus sequences carrying each allele '
            'observed in this run.'),
        'heterozygosity': (
            'The number and percentage of samples that are homozygous and '
            'heterozygous.')
    },
    name='Call alleles and genotypes',
    description=(
        'Call alleles and genotypes from per-sample consensus sequences '
        'using the Real Allele algorithm. Sequences below `min_coverage` '
        'reads are discarded, each sample keeps at most two sequences '
        '(subject to `proportion`), and identical sequences are grouped '
        'into alleles numbered by decreasing frequency (A1 is the most '
        'common). Provide `reference_alleles` to name alleles consistently '
        'across runs. Read counts come from the sequence IDs, or from '
        '`table` when it is provided.'),
    citations=[citations['Stone2024']],
    examples={'call_alleles_basic': call_alleles_basic,
              'call_alleles_with_reference': call_alleles_with_reference}
)
