************
Installation
************

``vcfsim`` is available for Linux and macOS. The instructions below cover
the standard conda install, plus an editable install from source for
development.

Via conda
=========

We recommend installing into a fresh conda environment. If you don't already
have a conda distribution installed, we recommend `Miniforge
<https://github.com/conda-forge/miniforge>`_ (lightweight, fully free,
conda-forge pre-configured, ships with the fast ``mamba`` solver) or
`Miniconda <https://docs.anaconda.com/miniconda/>`_::

    conda create -n vcfsim_env python=3.10
    conda activate vcfsim_env
    conda install -c bioconda vcfsim

Test the installation by running::

    vcfsim --help

For the full bioconda recipe and release history see
https://bioconda.github.io/recipes/vcfsim/README.html.

.. note::
    If the conda solver hangs on environment creation, install ``mamba`` and
    use it in place of ``conda`` for the install steps (``mamba install ...``).
    Miniforge ships with ``mamba`` already, and recent conda releases support
    ``--solver=libmamba`` for the same speedup.

From source
===========

To work from the latest code (or to contribute), clone the repository and
install in editable mode with ``pip``::

    git clone https://github.com/samuk-lab/vcfsim.git
    cd vcfsim
    pip install -e .
    vcfsim --help

This pulls in the runtime dependencies — ``numpy``, ``msprime``, and
``tskit`` — automatically. We recommend doing this inside a fresh conda
environment to keep the ``msprime`` build dependencies isolated from your
system Python.
