# ----------------------------------------------------------------------------
# Copyright (c) 2026, Oddant1.
#
# No license has been chosen for this software yet.
#
# See the LICENSE file distributed with this software.
# ----------------------------------------------------------------------------

import hashlib

import pandas as pd
import pandas.testing as pdt

import rachis
from rachis.plugin.testing import TestPluginBase
from rachis.plugin.util import transform
from q2_types.feature_data import DNAFASTAFormat

from q2_real_allele._methods import call_alleles


def _as_table(consensus):
    """Split consensus sequences into a feature table and its features.

    Features are named by the MD5 hash of their sequence, as dereplicated
    sequences usually are, so their IDs carry no sample or read count.
    """
    counts, features = {}, {}
    for id_, sequence in consensus.items():
        if 'CONSENSUS' in id_:
            continue
        sample, _, read_count = id_.rpartition('_')
        feature = hashlib.md5(str(sequence).encode()).hexdigest()
        features[feature] = str(sequence)
        counts.setdefault(sample, {})[feature] = int(read_count)
    # Samples stay in the order they first appear, because ties between
    # equally frequent alleles are named in the order they are first seen.
    table = pd.DataFrame(list(counts.values()), index=list(counts)).fillna(0)
    return table, pd.Series(features)


class RealAlleleTestCase(TestPluginBase):
    package = 'q2_real_allele.tests'

    def load_sequences(self, filename):
        return transform(self.get_data_path(filename),
                         from_type=DNAFASTAFormat, to_type=pd.Series)

    def import_sequences(self, filename):
        return rachis.Artifact.import_data('FeatureData[Sequence]',
                                           self.get_data_path(filename))

    def assertSequencesEqual(self, observed, expected):
        self.assertEqual([(id_, str(seq)) for id_, seq in observed.items()],
                         [(id_, str(seq)) for id_, seq in expected.items()])

    def assertMetadataEqual(self, observed, expected_filename):
        expected = rachis.Metadata.load(self.get_data_path(expected_filename))
        pdt.assert_frame_equal(observed.to_dataframe(),
                               expected.to_dataframe())

    def assertOutputsMatch(self, observed, expected_dir):
        alleles, genotyped, genotypes, frequencies, zygosity = observed
        self.assertSequencesEqual(
            alleles, self.load_sequences(expected_dir + '/alleles.fasta'))
        self.assertSequencesEqual(
            genotyped,
            self.load_sequences(expected_dir + '/genotyped-sequences.fasta'))
        self.assertMetadataEqual(genotypes,
                                 expected_dir + '/genotype-table.tsv')
        self.assertMetadataEqual(frequencies,
                                 expected_dir + '/allele-frequencies.tsv')
        self.assertMetadataEqual(zygosity,
                                 expected_dir + '/heterozygosity.tsv')


class CallAllelesTests(RealAlleleTestCase):
    def test_defaults(self):
        observed = call_alleles(self.load_sequences('consensus-1.fasta'))

        self.assertOutputsMatch(observed, 'expected/consensus-1')

    def test_reference_alleles_from_previous_run(self):
        first_alleles, *_ = call_alleles(
            self.load_sequences('consensus-1.fasta'))

        observed = call_alleles(self.load_sequences('consensus-2.fasta'),
                                reference_alleles=first_alleles)

        self.assertOutputsMatch(observed,
                                'expected/consensus-2-with-reference')

    def test_new_alleles_numbered_after_highest_reference_allele(self):
        alleles, _, _, frequencies, _ = call_alleles(
            self.load_sequences('consensus-1.fasta'),
            reference_alleles=self.load_sequences('reference-alleles.fasta'))

        # Reference alleles come first, in reference order, then new ones.
        self.assertEqual(list(alleles.index),
                         ['A7', 'A2', 'A8', 'A9', 'A10'])
        self.assertEqual(list(frequencies.to_dataframe().index),
                         ['A8', 'A2', 'A9', 'A10', 'A7'])

    def test_min_coverage(self):
        _, genotyped, genotypes, _, _ = call_alleles(
            self.load_sequences('consensus-1.fasta'), min_coverage=50)

        genotypes = genotypes.to_dataframe()
        # sample-05 has no sequence with 50 reads, and sample-08 keeps only
        # its 50-read sequence.
        self.assertNotIn('sample-05_S5', genotypes.index)
        self.assertEqual(list(genotypes.loc['sample-08_S8']),
                         [50, 'A3', 50, 'A3'])
        self.assertEqual(
            [id_ for id_ in genotyped.index if id_.startswith('sample-08')],
            ['sample-08_S8_50_A3'])

    def test_proportion_threshold_is_inclusive(self):
        # sample-08's second sequence has exactly 0.9 of the first's reads.
        *_, zygosity = call_alleles(
            self.load_sequences('consensus-1.fasta'), proportion=0.9)

        self.assertEqual(
            list(zygosity.to_dataframe()['sample-count']), [6, 1])

    def test_proportion_zero_keeps_two_sequences_per_sample(self):
        _, _, genotypes, _, zygosity = call_alleles(
            self.load_sequences('consensus-1.fasta'), proportion=0)

        self.assertEqual(
            list(genotypes.to_dataframe().loc['sample-03_S3']),
            [300, 'A1', 90, 'A3'])
        self.assertEqual(
            list(zygosity.to_dataframe()['sample-count']), [3, 4])

    def test_malformed_ids(self):
        with self.assertRaisesRegex(
                ValueError, r"<sample>_<read count>.*'sample-02_S2_many', "
                            r"'sample-03_S3_', '150'"):
            call_alleles(self.load_sequences('malformed-ids.fasta'))

    def test_no_sequences_pass_min_coverage(self):
        with self.assertRaisesRegex(
                ValueError, 'None of the 14 consensus sequences has at '
                            'least 1000 reads'):
            call_alleles(self.load_sequences('consensus-1.fasta'),
                         min_coverage=1000)

    def test_reference_allele_ids_must_be_numbered(self):
        with self.assertRaisesRegex(
                ValueError, r"A<number>.*'allele-2'"):
            call_alleles(
                self.load_sequences('consensus-1.fasta'),
                reference_alleles=self.load_sequences(
                    'reference-bad-ids.fasta'))

    def test_reference_allele_sequences_must_be_unique(self):
        with self.assertRaisesRegex(
                ValueError, 'share a sequence: A1 and A3'):
            call_alleles(
                self.load_sequences('consensus-1.fasta'),
                reference_alleles=self.load_sequences(
                    'reference-duplicate-sequences.fasta'))


class CallAllelesFromTableTests(RealAlleleTestCase):
    def test_matches_read_counts_from_ids(self):
        table, sequences = _as_table(
            self.load_sequences('consensus-1.fasta'))

        alleles, genotyped, genotypes, frequencies, zygosity = call_alleles(
            sequences, table=table)

        expected = 'expected/consensus-1/'
        self.assertSequencesEqual(
            alleles, self.load_sequences(expected + 'alleles.fasta'))
        self.assertMetadataEqual(frequencies,
                                 expected + 'allele-frequencies.tsv')
        self.assertMetadataEqual(zygosity, expected + 'heterozygosity.tsv')
        # The table is read one sample at a time rather than in the order of
        # the FASTA file, so only the order of these two differs.
        self.assertEqual(
            sorted((id_, str(seq)) for id_, seq in genotyped.items()),
            sorted((id_, str(seq)) for id_, seq in self.load_sequences(
                expected + 'genotyped-sequences.fasta').items()))
        pdt.assert_frame_equal(
            genotypes.to_dataframe().sort_index(),
            rachis.Metadata.load(self.get_data_path(
                expected + 'genotype-table.tsv')).to_dataframe().sort_index())

    def test_sequences_in_a_sample_can_share_a_read_count(self):
        sequences = pd.Series(
            {'f1': 'ACGTACGT', 'f2': 'TTGTACGT', 'f3': 'ACGAACGT'})
        table = pd.DataFrame(
            {'f1': [500, 1000], 'f2': [500, 500], 'f3': [0, 500]},
            index=['sample-a', 'sample-b'])

        alleles, genotyped, genotypes, _, _ = call_alleles(sequences,
                                                           table=table)

        self.assertEqual({id_: str(seq) for id_, seq in alleles.items()},
                         {'A1': 'ACGTACGT', 'A2': 'TTGTACGT'})
        # sample-b keeps f2 as its runner-up because it comes first in the
        # table, and f3 is dropped even though its ID would also have been
        # sample-b_500.
        self.assertEqual({id_: str(seq) for id_, seq in genotyped.items()}, {
            'sample-a_500_A1': 'ACGTACGT',
            'sample-a_500_A2': 'TTGTACGT',
            'sample-b_1000_A1': 'ACGTACGT',
            'sample-b_500_A2': 'TTGTACGT',
        })
        self.assertEqual(genotypes.to_dataframe().values.tolist(),
                         [[500, 'A1', 500, 'A2'], [1000, 'A1', 500, 'A2']])

    def test_table_features_must_be_sequences(self):
        with self.assertRaisesRegex(
                ValueError, "must be one of the sequences.*'f2'"):
            call_alleles(pd.Series({'f1': 'ACGTACGT'}),
                         table=pd.DataFrame({'f1': [20], 'f2': [30]},
                                            index=['sample-a']))


class RealAlleleV11RegressionTests(RealAlleleTestCase):
    """Compare with Real_AlleleV11 run on TEST_V1.fasta and TEST_V2.fasta.

    The input fixtures omit the repeated records of duplicated IDs (all with
    one or two reads), which QIIME 2 cannot import; Real_AlleleV11 produces
    identical output with and without them.
    """

    def read_v11_table(self, filename, **kwargs):
        return pd.read_csv(self.get_data_path('real-allele-v11/' + filename),
                           sep='\t', **kwargs)

    def test_matches_real_allele_v11(self):
        for name in 'TEST_V1', 'TEST_V2':
            with self.subTest(name):
                alleles, genotyped, genotypes, frequencies, zygosity = \
                    call_alleles(self.load_sequences(
                        'real-allele-v11/%s.fasta' % name))

                exp_alleles = self.read_v11_table(
                    name + '_alleles.tsv', header=None)
                self.assertEqual(
                    [(str(seq), id_) for id_, seq in alleles.items()],
                    list(exp_alleles.itertuples(index=False, name=None)))

                exp_frequencies = self.read_v11_table(
                    name + '_frequency.tsv', header=None)
                frequencies = frequencies.to_dataframe()['frequency']
                self.assertEqual(
                    [(str(alleles[id_]), count)
                     for id_, count in frequencies.items()],
                    list(exp_frequencies.itertuples(index=False, name=None)))

                self.assertSequencesEqual(
                    genotyped, self.load_sequences(
                        'real-allele-v11/%s.genotyped.fasta' % name))

                exp_genotypes = self.read_v11_table(
                    name + '_genotype_table.tsv', index_col=0)
                genotypes = genotypes.to_dataframe()
                self.assertEqual(list(genotypes.index),
                                 list(exp_genotypes.index))
                self.assertEqual(genotypes.values.tolist(),
                                 exp_genotypes.values.tolist())

                # Real_AlleleV11 miscounts zygosity, so derive the expected
                # counts from its genotype table instead of its report.
                heterozygous = int(
                    (exp_genotypes['Assigned allele'] !=
                     exp_genotypes['Assigned second allele']).sum())
                self.assertEqual(
                    list(zygosity.to_dataframe()['sample-count']),
                    [len(exp_genotypes) - heterozygous, heterozygous])


class CallAllelesActionTests(RealAlleleTestCase):
    def test_action_with_reference_from_previous_run(self):
        call = self.plugin.methods['call_alleles']

        first = call(sequences=self.import_sequences('consensus-1.fasta'))
        second = call(sequences=self.import_sequences('consensus-2.fasta'),
                      reference_alleles=first.alleles)

        for output in 'alleles', 'genotyped_sequences':
            self.assertEqual(str(getattr(second, output).type),
                             'FeatureData[Sequence]')
        for output in ('genotype_table', 'allele_frequencies',
                       'heterozygosity'):
            self.assertEqual(str(getattr(second, output).type),
                             'ImmutableMetadata')
        self.assertMetadataEqual(
            second.genotype_table.view(rachis.Metadata),
            'expected/consensus-2-with-reference/genotype-table.tsv')

    def test_action_with_table(self):
        call = self.plugin.methods['call_alleles']
        table, sequences = _as_table(
            self.load_sequences('consensus-1.fasta'))

        results = call(
            sequences=rachis.Artifact.import_data('FeatureData[Sequence]',
                                                  sequences),
            table=rachis.Artifact.import_data('FeatureTable[Frequency]',
                                              table))

        self.assertMetadataEqual(
            results.allele_frequencies.view(rachis.Metadata),
            'expected/consensus-1/allele-frequencies.tsv')

    def test_parameters_out_of_range(self):
        call = self.plugin.methods['call_alleles']
        sequences = self.import_sequences('consensus-1.fasta')

        for parameters in {'min_coverage': 0}, {'proportion': 1.5}:
            with self.subTest(**parameters):
                with self.assertRaisesRegex(TypeError, 'incompatible'):
                    call(sequences=sequences, **parameters)

    def test_examples(self):
        self.execute_examples()
