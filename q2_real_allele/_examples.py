# ----------------------------------------------------------------------------
# Copyright (c) 2026, Oddant1.
#
# No license has been chosen for this software yet.
#
# See the LICENSE file distributed with this software.
# ----------------------------------------------------------------------------

import importlib.resources

import rachis


def _consensus_sequences(filename):
    path = importlib.resources.files('q2_real_allele.tests') / 'data'
    return rachis.Artifact.import_data(
        'FeatureData[Sequence]', str(path / filename))


def consensus_sequences_1_factory():
    return _consensus_sequences('consensus-1.fasta')


def consensus_sequences_2_factory():
    return _consensus_sequences('consensus-2.fasta')


def _call_alleles(use, sequences, suffix='', **inputs):
    return use.action(
        use.UsageAction(plugin_id='real_allele', action_id='call_alleles'),
        use.UsageInputs(sequences=sequences, **inputs),
        use.UsageOutputNames(
            alleles='alleles' + suffix,
            genotyped_sequences='genotyped_sequences' + suffix,
            genotype_table='genotype_table' + suffix,
            allele_frequencies='allele_frequencies' + suffix,
            heterozygosity='heterozygosity' + suffix))


def call_alleles_basic(use):
    sequences = use.init_artifact('consensus_sequences',
                                  consensus_sequences_1_factory)

    alleles, _, genotype_table, *_ = _call_alleles(use, sequences)

    alleles.assert_output_type('FeatureData[Sequence]')
    genotype_table.assert_output_type('ImmutableMetadata')


def call_alleles_with_reference(use):
    sequences_1 = use.init_artifact('consensus_sequences_1',
                                    consensus_sequences_1_factory)
    sequences_2 = use.init_artifact('consensus_sequences_2',
                                    consensus_sequences_2_factory)

    use.comment('Call alleles in a first batch of samples.')
    alleles_1, *_ = _call_alleles(use, sequences_1, suffix='_1')

    use.comment('Name the alleles in a second batch consistently with the '
                'first by passing the first run\'s alleles as the reference. '
                'The new alleles output contains the reference alleles plus '
                'any new ones, ready to be used as the reference for the '
                'next batch.')
    alleles_2, *_ = _call_alleles(use, sequences_2, suffix='_2',
                                  reference_alleles=alleles_1)

    alleles_2.assert_output_type('FeatureData[Sequence]')
