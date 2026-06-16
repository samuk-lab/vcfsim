*********
Changelog
*********

Explanations of major changes to ``vcfsim`` are listed below. For
up-to-date info on minor versions and bugfixes, see the release notes on
GitHub: https://github.com/samuk-lab/vcfsim/releases

Unreleased
==========

Documentation
-------------

- Added Sphinx documentation hosted on Read the Docs.

vcfsim 1.0.29
=============

New features
------------

- **Comma-separated samples.** ``--samples`` and ``--samples_file`` now
  accept either comma- or whitespace-separated names, so the same list
  can be pasted into either flag.
- **Two-population mode** (``--population_mode 2``) for simulating a
  clean *C → A, B* split at a user-specified divergence time
  (``--div_time``).
- **HMM-based site missingness.** Four new flags (``--hmm_baseline``,
  ``--hmm_multiplier``, ``--hmm_p_low_to_high``, ``--hmm_p_high_to_low``)
  produce spatially clustered missingness instead of the uniform
  ``--percent_missing_sites`` pattern.

Bug fixes
---------

- Several VCF formatting fixes for downstream compatibility.

Packaging
---------

- Refactored VCF generation; bumped version and dependencies.
