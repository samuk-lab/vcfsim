************
About vcfsim
************

| **Maintainer:** Kieran Samuk (UC Riverside).
| **Original authors:** Paimon Goulart (UC Riverside) and Kieran Samuk (UC Riverside).

``vcfsim`` is a command-line tool that uses coalescent simulation to produce
realistic VCFs for use in methods development, benchmarking, and teaching. It
wraps `msprime <https://tskit.dev/msprime/>`_ for the underlying simulation
and adds a thin postprocessing layer that:

* Emits a single, well-formed VCF (with an ``msprime``-style header).
* Inserts **invariant sites** so the result is a true all-sites VCF.
* Applies configurable **missing data** at both the site and genotype level
  — either as uniform deterministic missingness, or as spatially clustered
  missingness governed by a two-state Hidden Markov Model.
* Supports arbitrary **ploidy**, custom **sample names**, and a simple
  two-population **demographic split**.
* Concatenates the output of multiple per-chromosome runs into one VCF
  when given a parameter file.

Why simulated all-sites VCFs?
=============================

Population genetic summary statistics — π, d\ :sub:`xy`, F\ :sub:`ST`,
Watterson's θ, Tajima's *D* — are sensitive to missing data, and many tools
silently treat missing genotypes as homozygous reference. The result is
biased estimates that are hard to spot without ground truth.

``vcfsim`` exists to provide that ground truth. Because each output VCF
comes from a coalescent simulation with known parameters (Ne, μ, ploidy,
divergence time, missingness rate), the *true* value of each summary
statistic is known by construction. That makes ``vcfsim`` a natural
companion for:

* **Benchmarking statistical estimators** — verify that an estimator
  recovers the correct value, and quantify the bias introduced by missing
  data.
* **Building test fixtures** — small, deterministic VCFs for unit and
  integration tests of downstream tools.
* **Teaching** — illustrate how missing data biases estimates of nucleotide
  diversity and divergence.
* **Methods development** — generate VCFs with controlled, spatially
  structured missingness to test mask-aware or block-aware estimators.

Notable features
================

* **All-sites VCF output** by default — both variant and invariant sites
  are written, so downstream tools that depend on callable-site denominators
  (e.g. `pixy <https://pixy.readthedocs.io>`_) work out of the box.
* **Two missing-data models.** A simple uniform model (``--percent_missing_sites``)
  and a two-state HMM that produces spatially clustered missingness
  (``--hmm_baseline``, ``--hmm_multiplier``, ``--hmm_p_low_to_high``,
  ``--hmm_p_high_to_low``).
* **Independent site- and genotype-level missingness** — control the rate
  at which whole rows go missing separately from the rate at which
  individual genotypes go missing within a row.
* **Arbitrary ploidy** via ``--ploidy``.
* **Custom sample names** from the command line (``--samples A1 B1 C1``)
  or from a file (``--samples_file names.txt``), with comma- or
  whitespace-separated input. Numeric ``--sample_size`` is the default.
* **Two-population mode** (``--population_mode 2``) that simulates an
  ancestral population *C* splitting into populations *A* and *B* at a
  user-specified divergence time (``--div_time``).
* **Multi-chromosome batch mode** — point ``--chromosome_file`` at a
  parameter file listing per-chromosome ploidy, sequence length, Ne, and
  mutation rate, and ``vcfsim`` runs each row and concatenates the
  results into a single VCF.
* **Reproducibility built in** — every run is keyed off a user-supplied
  ``--seed``. When ``--replicates > 1`` is requested, replicates are
  produced by incrementing the seed.

Acknowledgements
================

``vcfsim`` is built on top of `msprime <https://tskit.dev/msprime/>`_ and
`tskit <https://tskit.dev/>`_. The all-sites postprocessing and missing-data
models are designed to interoperate cleanly with the
`pixy <https://pixy.readthedocs.io>`_ population genetics toolkit.
