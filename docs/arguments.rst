*********
Arguments
*********

Below is the full list of arguments accepted by ``vcfsim``. You can also see
the canonical help text at any time with ``vcfsim --help``.

Required
========

**--seed [integer]**
    Random seed used to drive the coalescent simulation. The same seed (with
    the same parameters) reproduces the same VCF bit-for-bit. When
    ``--replicates > 1``, the seed is incremented by one for each replicate.

**--percent_missing_genotypes [integer]**
    Percent of *individual genotypes* (sample × site cells) to mark as
    missing in the VCF. Applied uniformly at random across the genotype
    matrix.

In addition, one of (sample specification)
==========================================

Exactly one of the following must be provided to specify the samples in the
VCF:

**--sample_size [integer]**
    Number of samples to draw from the simulated population. Sample IDs
    default to ``tsk_0``, ``tsk_1``, ..., ``tsk_n``.

**--samples [NAMES ...]**
    Custom sample names, comma- or space-separated (e.g. ``A1 B1 C1`` or
    ``A1,B1,C1``). The sample size is set to the number of names provided.

**--samples_file [path/to/names.txt]**
    Path to a file containing comma- or whitespace-separated sample names.

In addition, one of (missing-data model)
========================================

Exactly one of the following missing-data models must be specified. The
first produces uniform missingness across sites; the second produces
spatially clustered missingness from a two-state HMM.

**--percent_missing_sites [integer]**
    Percent of *sites* (VCF rows) to drop entirely. Missingness is applied
    uniformly at random across sites.

**HMM-based site missingness** *(all four flags required together)*
    These flags switch on a two-state Hidden Markov Model in which sites
    transition between a low-missing ("good") state and a high-missing
    ("bad") state. Sites in the bad state are missing at an elevated rate,
    producing spatially clustered patterns of missingness rather than
    uniform noise. The realized missingness fraction is stochastic but is
    approximately determined by the parameters.

    **--hmm_baseline [float, 0–1]**
        Baseline probability of site missingness when in the low-missing
        state.
    **--hmm_multiplier [float, ≥1]**
        Multiplier applied to ``--hmm_baseline`` when in the high-missing
        state.
    **--hmm_p_low_to_high [float, 0–1]**
        Per-site transition probability from low-missing to high-missing.
    **--hmm_p_high_to_low [float, 0–1]**
        Per-site transition probability from high-missing to low-missing.

Optional (simulation parameters)
================================

These describe the simulated chromosome and population. All are required
when running a single chromosome (i.e. when ``--chromosome_file`` is not
used), and all must be **omitted** when using ``--chromosome_file`` (the
per-chromosome values come from the file).

**--chromosome [string]**
    Chromosome name/label written into the ``#CHROM`` column of the VCF.

**--replicates [integer]**
    Number of replicate VCFs to produce. The seed is incremented by one for
    each replicate, so ``--seed 1234 --replicates 3`` produces three VCFs
    seeded with 1234, 1235, and 1236.

**--sequence_length [integer]**
    Length of the simulated chromosome in base pairs.

**--ploidy [integer]**
    Ploidy of the simulated samples (2 = diploid, 1 = haploid, etc.).

**--Ne [integer]**
    Effective population size of the simulated population.

**--mu [float, < 1]**
    Per-site, per-generation mutation rate.

Optional (output)
=================

**--output_file [prefix]**
    Prefix for the output VCF. The seed is appended automatically and the
    ``.vcf`` extension is added, so ``--output_file myvcf --seed 1234``
    writes to ``myvcf1234.vcf``. If omitted, the VCF is written to
    ``stdout``.

Optional (multi-chromosome batch mode)
======================================

**--chromosome_file [path/to/params.txt]**
    Path to a whitespace-separated file with one row per chromosome and
    five columns: ``chromosome ploidy sequence_length Ne mu``. Each row
    triggers a separate simulation and the results are concatenated into a
    single VCF. The numeric columns may be plain integers/floats or simple
    Python expressions (e.g. ``1e4``, ``1.7e6``, ``5.5e-9``). When
    ``--chromosome_file`` is used, the per-chromosome ``--chromosome``,
    ``--replicates``, ``--sequence_length``, ``--ploidy``, ``--Ne``, and
    ``--mu`` flags must all be omitted.

Optional (population structure)
===============================

**--population_mode [1|2]**
    Number of populations to simulate. ``1`` (default) simulates a single
    panmictic population. ``2`` simulates a split: an ancestral population
    *C* splits into *A* and *B* at ``--div_time`` generations before
    present. In mode 2, the number of samples must be even (they are split
    evenly between A and B).

**--div_time [integer]**
    Divergence time in generations before present for the *A*/*B* split.
    Required when ``--population_mode 2``; ignored otherwise.

Help
====

**--help**
    Print the full help message and exit.

Example
=======

A typical single-chromosome run with 10 diploid samples, no missing data,
and ten thousand simulated sites:

.. code:: console

    vcfsim \
        --chromosome 1 \
        --replicates 1 \
        --seed 1234 \
        --sequence_length 10000 \
        --ploidy 2 \
        --Ne 100000 \
        --mu 1e-6 \
        --percent_missing_sites 0 \
        --percent_missing_genotypes 0 \
        --sample_size 10 \
        --output_file myvcf

This writes the result to ``myvcf1234.vcf``.
