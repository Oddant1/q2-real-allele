# q2-real-allele

A [QIIME 2](https://qiime2.org) plugin that calls alleles and genotypes from per-sample consensus amplicon sequences, using the [Real Allele](https://github.com/GrantPem/Real_Allele) algorithm by Jason Sahl and Grant Pemberton (packaged from `Real_AlleleV11`).
Every file that `Real_AlleleV11` reads or writes is replaced by a QIIME 2 artifact of an existing type from `q2-types`.

## Installation instructions

**The following instructions are intended to be a starting point** and should be replaced when `q2-real-allele` is ready to share with others.
They will enable you to install the most recent *development* version of `q2-real-allele`.
Remember that *release* versions should be used for all "real" work (i.e., where you're not testing or prototyping) - if there aren't instructions for installing a release version of this plugin, it is probably not yet intended for use in practice.

### Install Prerequisites

[Miniconda](https://conda.io/miniconda.html) provides the `conda` environment and package manager, and is currently the only supported way to install QIIME 2.
Follow the instructions for downloading and installing Miniconda.

After installing Miniconda and opening a new terminal, make sure you're running the latest version of `conda`:

```bash
conda update conda
```

###  Install development version of `q2-real-allele`

Next, you need to get into the top-level `q2-real-allele` directory.
If you already have this (e.g., because you just created the plugin), this may be as simple as running `cd q2-real-allele`.
If not, you'll need the `q2-real-allele` directory on your computer.
How you do that will differ based on how the package is shared, and ideally the developer will update these instructions to be more specific (remember, these instructions are intended to be a starting point).
For example, if it's maintained in a GitHub repository, you can achieve this by [cloning the repository](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository).
Once you have the directory on your computer, change (`cd`) into it.

If you're in a conda environment, deactivate it by running `conda deactivate`.


Then, follow the install instructions below, based on your machine's architecture:

<details>
<summary><strong>🍏&nbsp;Apple Silicon (ARM)</strong></summary>
<p>&nbsp;</p>

Start by creating a new conda environment:

```shell
CONDA_SUBDIR=osx-64 conda env create -n q2-real-allele-dev --file ./environment-files/q2-real-allele-rachis-tiny-dev.yml
```

After this completes, activate the new environment you created by running:

```shell
conda activate q2-real-allele-dev
```

Once this new environment has been activated, update your conda config to set the subdir to osx-64:

```shell
conda config --env --set subdir osx-64
```

Finally, run:

```shell
make install
```
</details>

<details>
<summary><strong>🛠&nbsp;All other architectures (Apple Intel, Linux, WSL)</strong></summary>
<p>&nbsp;</p>

Start by creating a new conda environment:

```shell
conda env create -n q2-real-allele-dev --file ./environment-files/q2-real-allele-rachis-tiny-dev.yml
```

After this completes, activate the new environment you created by running:

```shell
conda activate q2-real-allele-dev
```

Finally, run:

```shell
make install
```
</details>

## Testing the most recent development version of `q2-real-allele`

After completing the install steps above, confirm that everything is working as expected by running:

```shell
make test
```

You should get a report that tests were run, and you should see that all tests passed and none failed.
It's usually ok if some warnings are reported.

If all of the tests pass, make QIIME 2's command line interface aware of `q2-real-allele` by running:

```shell
qiime dev refresh-cache
```

You should then see the plugin in the list of available plugins if you run `qiime info`, and you can review the help text by running:

```shell
qiime real-allele call-alleles --help
```

## Usage

### 1. Import your consensus sequences

Each sequence ID must have the form `<sample>_<read count>`: everything before the last underscore identifies the sample, and the final field is the number of reads supporting that consensus sequence.
For example, `RHIMI-ARA13-xx-DN-020-JB_S86_17441` is sample `RHIMI-ARA13-xx-DN-020-JB_S86` with 17,441 reads.
Sequences whose ID contains `CONSENSUS` are ignored.

QIIME 2 requires sequence IDs to be unique.
Consensus FASTA files can contain duplicate IDs when two low-abundance sequences from the same sample have the same read count (e.g. `TEST_V1.fasta` contains 67 of these, all with 1 or 2 reads).
Because sequences below `--p-min-coverage` reads are discarded anyway, you can remove them before importing without changing the results.
Use the same threshold for both steps (10 is the default):

```shell
awk -v min=10 '/^>/ {n = split($1, f, "_"); keep = (f[n] + 0 >= min)} keep' \
  TEST_V1.fasta > TEST_V1.filtered.fasta

qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path TEST_V1.filtered.fasta \
  --output-path TEST_V1.qza
```

### 2. Call alleles

```shell
qiime real-allele call-alleles \
  --i-sequences TEST_V1.qza \
  --output-dir TEST_V1-results
```

`--p-min-coverage` (default 10) and `--p-proportion` (default 0.5) are the `-c` and `-p` options of `Real_AlleleV11`.

### Keeping rare alleles

`Real_AlleleV11` keeps at most two sequences per sample: the most abundant, and the second most abundant if it has at least `--p-proportion` as many reads.
Pass `--p-keep-rare-alleles` to call alleles as `Rare_Allele` does instead, keeping every sequence with at least `--p-proportion` as many reads as the sample's most abundant, so that a sample can have any number of alleles:

```shell
qiime real-allele call-alleles \
  --i-sequences TEST_V1.qza \
  --p-keep-rare-alleles \
  --output-dir TEST_V1-rare-results
```

A sample's two most abundant alleles are the ones `Real_AlleleV11` would have called, and any after them are its rare alleles.
The genotype table gets a read count and allele column for each, `read-count-3` and `allele-3` onwards, with as many as the sample with the most alleles needs.
As in `Rare_Allele`, each allele is listed once, so a homozygous sample has values only in `read-count-1` and `allele-1`, and a sample with fewer alleles than another has no value in the columns it lacks.
In the heterozygosity summary, a sample with more than one allele counts as heterozygous.

### Giving read counts in a feature table

Instead of encoding read counts in the sequence IDs, you can pass a `FeatureTable[Frequency]` with `--i-table`.
Its features are the `--i-sequences` (whose IDs can then be anything, such as MD5 hashes), its samples are the samples to genotype, and each value is the number of reads supporting that sequence in that sample.
Because read counts no longer have to be unique within a sample, two sequences in a sample can have the same count, such as the two alleles of a heterozygote.

The `dereplicated_tables` and `dereplicated_sequences` outputs of [q2-multi-amplicon](https://example.com/q2-multi-amplicon)'s `multi-amplicon-analysis` and `dereplicate-classifications` have this form, with one table and one set of sequences per amplicon:

```shell
qiime real-allele call-alleles \
  --i-sequences dereplicated-sequences/ref1.qza \
  --i-table dereplicated-tables/ref1.qza \
  --output-dir ref1-results
```

The genotyped sequences are still named `<sample>_<read count>_<allele>`.
Records are read one sample at a time in the order of the table, so where alleles are equally frequent, they are numbered in the order their samples appear in the table.

### 3. Name alleles consistently across runs

Pass the `alleles` output of an earlier run as `--i-reference-alleles`.
Sequences identical to a reference allele keep its name, and new alleles are numbered after the highest reference allele.
The new `alleles` output contains every reference allele plus the new ones, so it can be the reference for the next run.

```shell
qiime real-allele call-alleles \
  --i-sequences TEST_V2.qza \
  --i-reference-alleles TEST_V1-results/alleles.qza \
  --output-dir TEST_V2-results
```

To reuse an `_alleles.tsv` dictionary written by `Real_AlleleV11` or `allele_driverV5.py`, convert it to FASTA and import it (this also adds the `A` prefix to the numeric allele names written by `allele_driverV5.py`):

```shell
awk -F'\t' '{id = $2; if (id !~ /^A/) id = "A" id; print ">" id "\n" $1}' \
  TEST_V1_alleles.tsv > TEST_V1_alleles.fasta

qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path TEST_V1_alleles.fasta \
  --output-path TEST_V1-reference-alleles.qza
```

### 4. View the results

The tables are `ImmutableMetadata` artifacts, so you can view them with `qiime metadata tabulate` or export them as TSV files:

```shell
qiime metadata tabulate \
  --m-input-file TEST_V1-results/genotype_table.qza \
  --o-visualization genotype-table.qzv

qiime tools export \
  --input-path TEST_V1-results/genotype_table.qza \
  --output-path genotype-table
```

## Inputs and outputs

| `Real_AlleleV11` file | `call-alleles` artifact | QIIME 2 type | Contents |
| --- | --- | --- | --- |
| `-f <name>.fasta` | `--i-sequences` | `FeatureData[Sequence]` | Consensus sequences to genotype. |
| `-a <name>_alleles.tsv` | `--i-reference-alleles` | `FeatureData[Sequence]` | Previously named alleles (IDs `A1`, `A2`, ...). |
| — | `--i-table` | `FeatureTable[Frequency]` | Optional read counts of each sequence in each sample, in place of the counts in the sequence IDs. |
| `<name>_alleles.tsv` | `--o-alleles` | `FeatureData[Sequence]` | Allele sequences with their names as IDs. |
| `<name>.genotyped.fasta` | `--o-genotyped-sequences` | `FeatureData[Sequence]` | Retained consensus sequences, with the allele appended to the ID. |
| `<name>_genotype_table.tsv` | `--o-genotype-table` | `ImmutableMetadata` | Per sample: `read-count-1`, `allele-1`, `read-count-2`, `allele-2`. |
| `<name>_frequency.tsv` | `--o-allele-frequencies` | `ImmutableMetadata` | Per allele: the number of retained sequences carrying it. |
| `<name>_Heterozygosity.txt` | `--o-heterozygosity` | `ImmutableMetadata` | Number and percentage of homozygous and heterozygous samples. |

Tables that `Real_AlleleV11` keyed by sequence are keyed by allele name here.
The sequences themselves are in the `alleles` artifact.

## Differences from `Real_AlleleV11`

When no reference alleles are provided, the alleles, genotyped sequences, genotype table and allele frequencies match those written by `Real_AlleleV11` (this is checked by the test suite using `TEST_V1.fasta` and `TEST_V2.fasta`).
These behaviors were changed to fix bugs:

- **Reference alleles work.** `Real_AlleleV11 -a` stopped with a `NameError` because `load_previous_alleles` was never defined.
- **The alleles output keeps allele names.** With `-a`, `Real_AlleleV11` overwrote `<name>_alleles.tsv` with sequence counts, so it could not be used as the next reference.
- **The genotype table uses the assigned allele names.** With `-a`, `Real_AlleleV11` renumbered alleles by frequency in the genotype table, disagreeing with the genotyped FASTA.
- **Heterozygosity is counted correctly.** `Real_AlleleV11` counted every sample as a homozygote, and missed heterozygotes whose more abundant sequence came second in the file. For `TEST_V1.fasta` it reported 86 homozygotes and 1.15% heterozygosity instead of 85 and 1.16%. The table is also produced when there are no heterozygotes.
- **Malformed IDs are errors.** IDs without a numeric read count after the last underscore are reported instead of being skipped or crashing the run.

## About

The `q2-real-allele` Python package was [created from a template](https://develop.qiime2.org/en/stable/plugins/tutorials/create-from-template.html).
To learn more about `q2-real-allele`, refer to the [project website](https://example.com).
To learn how to use QIIME 2, refer to the [QIIME 2 User Documentation](https://docs.qiime2.org).
To learn QIIME 2 plugin development, refer to [*Developing with QIIME 2*](https://develop.qiime2.org).

If you use this plugin, please cite:

Stone NE, Ballard R, Bourgeois RM, Pemberton GL, et al. (2024) A mutation associated with resistance to synthetic pyrethroids is widespread in US populations of the tropical lineage of *Rhipicephalus sanguineus* s.l. *Ticks and Tick-borne Diseases* 15(4):102344. https://doi.org/10.1016/j.ttbdis.2024.102344

`q2-real-allele` is a QIIME 2 community plugin, meaning that it is not necessarily developed and maintained by the developers of QIIME 2.
Please be aware that because community plugins are developed by the QIIME 2 developer community, and not necessarily the QIIME 2 developers themselves, some may not be actively maintained or compatible with current release versions of the QIIME 2 distributions.
More information on development and support for community plugins can be found [here](https://library.qiime2.org).
If you need help with a community plugin, first refer to the [project website](https://example.com).
If that page doesn't provide information on how to get help, or you need additional help, head to the [Community Plugins category](https://forum.qiime2.org/c/community-contributions/community-plugins/14) on the QIIME 2 Forum where the QIIME 2 developers will do their best to help you.

## Standalone scripts

The original command-line scripts are still in the top-level directory and are not part of the plugin.
Their documentation follows unchanged.

### Real_Allele
-Uses biopython version 1.81

- The Real_Allele Script is designed to determine allele frequencies and assign allele numbers from a consensus FASTA file containing DNA sequences.  
  
  
- Sequencing artifacts will be filtered out based on a criterion of zygosity proportion, and minimum read count.
  
- Samples with more than two alleles will filter out the sequences with alleles lower than the top two read counts.

- Alleles are assigned based on their frequency, allele 1 will be assigned to the most common allele.  
- Creates a TSV dictionary file containing sequences and their counts.

  
- The FASTA file is needed to be named in this format: Gene_Exon.fasta (A minimum of 2 fields with underscores is required)
- The sequences within the file also need to be named to this convention:  
RHIMI-ARA13-xx-DN-020-JB_S86_17441  
"RHIMI-," "-xx-DN-020-JB", "-xx-VE-020-JB" will be removed.  
ARA13 is the sample name.  
17441 is the read count. 
This assumes three fields, all separated by "_," with the last field being the read count.  
The output name will have five elements: ARA13_AQP1_E0101_17441_A1.

- The Real_Allele will export additional CSV helpful for analysis.  
alleles.tsv = The found alleles and the assigned allele number.  
frequency.tsv = The total number of sequences for each allele number.  
genotype_table.tsv = A table showing each sequences assigned allele number, and total read count for each. (Real_Allele will always show two, if it is homozygote it will have allele 1 in both columns, Rare alleles will show all alleles that pass the filter not just the top two.)  
Heterozygosity.txt = Will show the total count of heterozygotes and homozygotes in the fasta file, and will printe the pergencatge of heterzygosity.   
  
  
- Test files are included with the repository. First, run the script on an initial file:  
python allele_driver.py -f TEST_V1.fasta  
--a command will use the TSV dictionary file and apply the previous allele naming scheme to the new fasta file, keeping the alleles defined in the dictionary uniform:  
python allele_driver.py -f TEST_V2.fasta -a TEST_V1_alleles.tsv  
Help command python allele_driverV5.py -h


### Rare_Allele 
-Uses biopython version 1.81

- The Rare_Allele Script is designed to determine allele frequencies and assign allele numbers from a consensus FASTA file containing DNA sequences while still retaining all sequences that pass the filter.

- Can have an indefinite number of possible alleles assigned to one sample  

### Real_Allele_Pipeline & Rare_Allele_Piepline
- The Real/Rare_Allele_Pipeline automatically extracts raw genomic sequences filters, transforms, and integrate the data for genotype analysis.
- The pipeline will first reoncile the dna sequence orientation, then it will prompt the user to remove the primer sequence from a fasta file. After that it will then run real_allele then filter out redundant seuqneces using allele_cleaner.
- Only has one manual step, the system will prompt the user to remove primers from the sequence.
- I will work on a method to remove primers automatically without user input. 

### Real_Allele 5
- allele_driverV5.py is the version of the script used in this publication: https://www.sciencedirect.com/science/article/pii/S1877959X24000372
 
