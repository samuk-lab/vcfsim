.. vcfsim documentation master file.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. raw:: html

    <div align="center"><h1> vcfsim 1.0 </h1></div>

.. image:: images/vcfsim_logo.png
   :width: 200
   :align: center

What is vcfsim?
===============
``vcfsim`` is a command-line tool for generating simulated VCFs (Variant Call Format files used to encode genetic variation data). It pairs a coalescent simulation backend (`msprime <https://tskit.dev/msprime/>`_) with lightweight postprocessing to produce biologically realistic VCFs with parameterized missing data — a full simulated dataset can be created from just a few command-line arguments.

In particular, ``vcfsim`` makes it easy to simulate **all-sites VCFs** that contain both variant and invariant sites. All-sites VCFs are required for unbiased estimation of π and d\ :sub:`xy` (see `pixy <https://pixy.readthedocs.io>`_), and are typically expensive to obtain from real data — ``vcfsim`` is designed to drop straight into pixy-style workflows for testing, benchmarking, and methods development.

``vcfsim`` also supports two missing-data models (uniform and HMM-based spatial clustering), arbitrary ploidy, custom sample names, two-population splits, and multi-chromosome batch runs from a parameter file.

.. toctree::
   :caption: Documentation
   :maxdepth: -1

   about
   installation
   arguments
   changelog
   contributing


.. toctree::
    :maxdepth: -1
    :caption: Guides

    usage
    missing_data
    populations
    multi_chromosome

How should I cite vcfsim?
=========================
If you use ``vcfsim`` in your research, please cite the repository and the underlying coalescent simulator:

* ``vcfsim``: https://github.com/samuk-lab/vcfsim
* Baumdicker, F. *et al.* (2022). Efficient ancestry and mutation simulation with msprime 1.0. *Genetics*, 220(3), iyab229. https://doi.org/10.1093/genetics/iyab229
