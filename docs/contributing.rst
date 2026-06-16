************
Contributing
************

Contributions to ``vcfsim`` are welcome — bug reports, documentation
fixes, new features, and improved tests are all useful. This page covers
how to set up a development environment, run the tool from source, and
open a pull request.

If you are just looking to install ``vcfsim`` for use rather than
development, see :doc:`installation`.

Installing the development version
==================================

Clone the repository and install in editable mode. We recommend doing this
inside a fresh conda environment so the ``msprime`` build dependencies
are isolated from your system Python::

    git clone https://github.com/samuk-lab/vcfsim.git
    cd vcfsim
    conda create -n vcfsim_dev python=3.10
    conda activate vcfsim_dev
    pip install -e .

You should now be able to run the development version directly::

    vcfsim --help

Project layout
==============

A quick map of the repository for orientation:

* ``vcfsim/`` — the main package.

  * ``__main__.py`` — CLI entry point and argument wiring.
  * ``SimulatorClass.py`` — ``msprime`` setup, postprocessing, and VCF
    emission.

* ``docs/`` — Sphinx documentation source (this page lives here).
* ``setup.py`` — package metadata and entry point.
* ``meta.yaml`` — bioconda recipe.

Working on the documentation
============================

The documentation lives under ``docs/`` and is built with Sphinx. To
preview your changes locally, install the doc dependencies and run a
build::

    pip install -r docs/requirements.txt
    sphinx-build -b html docs/ /tmp/vcfsim_html

Then open ``/tmp/vcfsim_html/index.html`` in a browser. On Windows,
substitute a writable path such as ``%TEMP%\\vcfsim_html`` and open the
resulting ``index.html`` with ``start``.

Reporting bugs and requesting features
======================================

Bug reports, feature requests, and questions should go to the
`GitHub issue tracker <https://github.com/samuk-lab/vcfsim/issues>`_.
When reporting a bug, the most useful information to include is:

* Your operating system and Python version.
* The full ``vcfsim`` command that triggered the problem.
* The full error output (please paste in text, not screenshots, when
  possible).
* The chromosome parameter file, if you were using ``--chromosome_file``.

Because every ``vcfsim`` run is keyed off a ``--seed``, including the
exact command is usually enough to reproduce the issue exactly — no need
to attach output VCFs.

Pull request workflow
=====================

If you're new to contributing on GitHub, the basic flow is:

1. Fork ``samuk-lab/vcfsim`` on GitHub.
2. Clone your fork locally and create a topic branch::

       git checkout -b fix/missing-data-edge-case

3. Make your changes and run a quick smoke test::

       vcfsim --seed 1 --chromosome 1 --replicates 1 \
              --sequence_length 1000 --ploidy 2 --Ne 10000 --mu 1e-6 \
              --percent_missing_sites 0 --percent_missing_genotypes 0 \
              --sample_size 4

4. Commit with a clear message describing *what* changed and *why*.
5. Push the branch to your fork and open a pull request against
   ``samuk-lab/vcfsim:main``.

A good pull request:

* **Has a descriptive title** that summarizes the change in plain English.
* **Links any related issues** in the description (``Fixes #123``).
* **Stays focused** on a single concern — one PR per logical change makes
  review much faster.
* **Updates the docs** under ``docs/`` if a user-facing flag or output
  changes.

License
=======

``vcfsim`` is released under the MIT License (see ``LICENSE.txt`` in the
repo root). By contributing code, documentation, or other materials to
the project, you agree that your contributions will be made available
under the same license.
