*******************
Missing data models
*******************

``vcfsim`` offers two independent levers for controlling missing data:

* **Genotype-level missingness** (``--percent_missing_genotypes``) — a
  per-cell drop rate applied uniformly across the genotype matrix. Always
  required, regardless of the site-missingness model.
* **Site-level missingness** — controlled by either a uniform model
  (``--percent_missing_sites``) or a Hidden Markov Model
  (``--hmm_*`` flags). Exactly one of these must be specified.

The two layers compose: site-level missingness drops entire VCF rows; the
genotype-level rate then applies to the cells of the rows that remain.

Uniform site missingness
========================

The simplest model drops a fixed percentage of VCF rows uniformly at
random:

.. code:: console

    vcfsim \
        --chromosome 1 --replicates 1 --seed 1234 \
        --sequence_length 10000 --ploidy 2 --Ne 100000 --mu 1e-6 \
        --percent_missing_sites 20 \
        --percent_missing_genotypes 5 \
        --sample_size 10 \
        --output_file myvcf

Here roughly 20% of sites are dropped, and within the surviving sites
roughly 5% of genotypes are marked missing. The realized percentages are
close to but not exactly the requested values for small sample sizes.

HMM-based spatially clustered missingness
=========================================

Real sequencing data rarely produces uniformly random missing sites —
poorly mappable regions and low-coverage stretches produce *runs* of
missing data instead. To simulate this pattern, leave
``--percent_missing_sites`` unset and supply all four HMM parameters:

.. code:: console

    vcfsim \
        --chromosome 1 --replicates 1 --seed 4000 \
        --sequence_length 10000 --ploidy 2 --Ne 100000 --mu 1e-6 \
        --percent_missing_genotypes 0 \
        --hmm_baseline 0.05 \
        --hmm_multiplier 6 \
        --hmm_p_low_to_high 0.002 \
        --hmm_p_high_to_low 0.005 \
        --sample_size 10 \
        --output_file myvcf

How the HMM works
-----------------

Each site is in one of two hidden states:

* **Low-missing** ("good") — the per-site missingness probability is
  ``--hmm_baseline``.
* **High-missing** ("bad") — the per-site missingness probability is
  ``--hmm_baseline × --hmm_multiplier``.

The state evolves along the chromosome as a two-state Markov chain with
transition probabilities ``--hmm_p_low_to_high`` (low → high) and
``--hmm_p_high_to_low`` (high → low). The expected fraction of sites in
the high-missing state is approximately
``p_low_to_high / (p_low_to_high + p_high_to_low)``; combined with the
two per-state missingness probabilities, this determines the long-run
overall missing-site rate.

The expected mean length of a high-missing run is approximately
``1 / p_high_to_low`` sites, so smaller ``--hmm_p_high_to_low`` values
produce longer contiguous bad regions.

Picking parameters
------------------

The example above produces clearly clustered missingness without
overwhelming the dataset:

* Baseline 5% missingness in good regions.
* 6× elevated rate (30%) in bad regions.
* Bad regions are entered on average once per ~500 sites
  (``1 / 0.002``).
* Bad regions persist on average ~200 sites (``1 / 0.005``).

To increase clustering without changing the overall missing rate, lower
both transition probabilities by the same factor. To increase the overall
missing rate while keeping the same spatial structure, raise
``--hmm_baseline`` and/or ``--hmm_multiplier``.

Mixing the two models
=====================

If both ``--percent_missing_sites`` and HMM flags are passed, ``vcfsim``
warns and uses the uniform model. To use the HMM, leave
``--percent_missing_sites`` unset.
