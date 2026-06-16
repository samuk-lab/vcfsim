******************************
Simulating a population split
******************************

By default ``vcfsim`` simulates a single panmictic population. To produce a
VCF that contains two populations with a shared history, switch to
two-population mode:

.. code:: console

    vcfsim \
        --chromosome 1 --replicates 1 --seed 1234 \
        --sequence_length 10000 --ploidy 2 --Ne 100000 --mu 1e-6 \
        --percent_missing_sites 0 --percent_missing_genotypes 0 \
        --sample_size 10 \
        --population_mode 2 --div_time 1000 \
        --output_file myvcf

The model
=========

``--population_mode 2`` simulates a clean two-population split: an
ancestral population *C* of effective size ``--Ne`` splits into two
present-day populations *A* and *B*, both of size ``--Ne``, at
``--div_time`` generations before present.

* **--population_mode 2** activates the split.
* **--div_time** sets the split time in generations before present.
  Required when ``--population_mode 2``; ignored otherwise.

Sample distribution
===================

Samples are split evenly between the two populations, so the total sample
count must be even. This applies to all three sample-specification flags:

* ``--sample_size 10`` puts 5 samples in *A* and 5 in *B*.
* ``--samples A1 A2 A3 A4`` puts the first two in *A* and the last two in
  *B*.
* ``--samples_file names.txt`` likewise splits the names in order.

If an odd number of samples is requested in mode 2, ``vcfsim`` raises an
error rather than silently truncating.

Using the output for F\ :sub:`ST` simulations
=============================================

Two-population mode is the natural choice for benchmarking
between-population statistics (F\ :sub:`ST`, d\ :sub:`xy`). The split time
gives you a single, intuitive knob for the level of differentiation:
larger ``--div_time`` produces deeper splits and higher F\ :sub:`ST`.

When using the output for `pixy <https://pixy.readthedocs.io>`_-style
analyses, build the populations file from the sample names emitted by
``vcfsim`` — the first half belong to population *A* and the second half
to population *B*.
