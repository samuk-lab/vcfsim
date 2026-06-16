***************
Basic usage
***************

This page walks through the most common invocations of ``vcfsim``. For the
full argument reference, see :doc:`arguments`.

A first run
===========

The minimal command needed to produce a VCF specifies a seed, a
missing-data model (here ``--percent_missing_sites 0`` and
``--percent_missing_genotypes 0`` for a clean reference dataset), and a
way to specify samples:

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
        --output_file myvcf \
        --sample_size 10

This produces ``myvcf1234.vcf`` — the seed is appended to the prefix
automatically. If ``--output_file`` is omitted, the VCF is written to
``stdout`` instead, which is convenient for piping into other tools.

Producing replicates
====================

Asking for more than one replicate runs the simulator multiple times with
an incrementing seed. With ``--seed 1234 --replicates 3``, ``vcfsim``
writes ``myvcf1234.vcf``, ``myvcf1235.vcf``, and ``myvcf1236.vcf``. This
is the recommended way to produce a set of independent replicates with
known seeds — it keeps the mapping from seed to file deterministic and
inspectable.

Custom sample names
===================

By default samples are named ``tsk_0``, ``tsk_1``, ..., ``tsk_n``. To
write explicit names into the VCF columns, pass ``--samples`` with the
names directly:

.. code:: console

    vcfsim \
        --chromosome 1 --replicates 1 --seed 1234 \
        --sequence_length 10000 --ploidy 2 --Ne 100000 --mu 1e-6 \
        --percent_missing_sites 0 --percent_missing_genotypes 0 \
        --samples A1 B1 C1 D1 \
        --output_file myvcf

Names may be space-separated (as above) or comma-separated
(``--samples A1,B1,C1,D1``). The sample size is set to the number of
names provided — ``--sample_size`` and ``--samples`` are mutually
exclusive.

For larger sample sets, store the names in a file and use
``--samples_file``:

.. code:: console

    vcfsim \
        --chromosome 1 --replicates 1 --seed 1234 \
        --sequence_length 10000 --ploidy 2 --Ne 100000 --mu 1e-6 \
        --percent_missing_sites 0 --percent_missing_genotypes 0 \
        --samples_file names.txt \
        --output_file myvcf

where ``names.txt`` contains comma- or whitespace-separated names::

    A1 B1 C1 D1 E1

Writing to stdout
=================

Omitting ``--output_file`` sends the VCF to ``stdout``, so you can pipe
the simulator directly into ``bgzip``, ``bcftools``, or any other VCF
consumer:

.. code:: console

    vcfsim \
        --chromosome 1 --replicates 1 --seed 1234 \
        --sequence_length 10000 --ploidy 2 --Ne 100000 --mu 1e-6 \
        --percent_missing_sites 0 --percent_missing_genotypes 0 \
        --sample_size 10 \
        | bgzip -c > myvcf.vcf.gz

This is convenient for one-off runs and is the form to use when feeding
``vcfsim`` output into a downstream pipeline without leaving intermediate
files on disk.
