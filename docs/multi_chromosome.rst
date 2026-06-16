**********************************
Multi-chromosome runs
**********************************

For genome-scale simulations — or any time you want to vary parameters
across chromosomes — pass a parameter file to ``--chromosome_file``
instead of specifying the per-chromosome parameters on the command line.
``vcfsim`` runs each row of the file as an independent simulation and
concatenates the results into a single VCF.

The parameter file format
=========================

The file is whitespace-separated, headerless, and has exactly five
columns::

    chromosome   ploidy   sequence_length   Ne   mu

For example::

    1   2   10000   100000   1e-6
    2   2   15000   100000   1e-6
    X   1    8000    50000   1e-6

* ``chromosome`` — name written into the ``#CHROM`` column for the
  corresponding rows.
* ``ploidy`` — must be a positive integer.
* ``sequence_length``, ``Ne``, ``mu`` — may be plain numbers or simple
  Python numeric expressions (``1e4``, ``1.7e6``, ``5.5e-9``, ``2*10000``).

Each row is run as a separate ``msprime`` simulation; the per-chromosome
ploidy, sequence length, Ne, and mutation rate can therefore vary freely.
This is particularly useful for simulating sex chromosomes alongside
autosomes — give the X-chromosome row ``ploidy = 1`` and the autosome
rows ``ploidy = 2``.

Running with a parameter file
=============================

When using ``--chromosome_file``, omit the per-chromosome flags
(``--chromosome``, ``--replicates``, ``--sequence_length``, ``--ploidy``,
``--Ne``, ``--mu``) — they come from the file instead. Everything else
(``--seed``, sample specification, missingness model, ``--output_file``)
works as usual:

.. code:: console

    vcfsim \
        --seed 1234 \
        --percent_missing_sites 0 \
        --percent_missing_genotypes 0 \
        --output_file myvcf \
        --sample_size 10 \
        --chromosome_file input.txt

This writes a single concatenated VCF (``myvcf.vcf``) containing all
chromosomes from ``input.txt``. The header is written once at the top
and the per-chromosome bodies are appended in order.

Combining with custom sample names
==================================

The sample-specification flags (``--sample_size``, ``--samples``,
``--samples_file``) work in batch mode exactly as they do for single-chromosome
runs. To label the same set of samples across all chromosomes:

.. code:: console

    vcfsim \
        --seed 1234 \
        --percent_missing_sites 0 \
        --percent_missing_genotypes 0 \
        --output_file myvcf \
        --samples_file names.txt \
        --chromosome_file input.txt

The same sample names are used in every per-chromosome simulation, so the
final concatenated VCF has a single, consistent column ordering.
