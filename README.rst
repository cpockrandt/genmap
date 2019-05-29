GenMap 1.0 |buildstatus|
------------------------

Fast and Exact Computation of Genome Mappability
================================================

.. |BUILDSTATUS| image:: https://travis-ci.org/cpockrandt/genmap.svg?branch=master
    :target: https://travis-ci.org/cpockrandt/genmap

GenMap is a tool to compute the mappability respectively frequency of nucleotide sequences (DNA and RNA).
In particular, it computes the (k,e)-frequency, i.e., how often each k-mer from the sequence occurs with up to e errors
in the sequence itself.
The (k,e)-mappability is the inverse of the (k,e)-frequency.
Hence, a mappability value of 1 at position i indicates that the k-mer in the sequence at position i occurs only once
in the sequence with up to e errors.
A low mappability value indicates that this k-mer belongs to a repetitive region.

A small example on how to run GenMap is listed below, for detailed examples such as marker sequence computation on
multiple fasta files, please check out our GitHub `Wiki pages <https://github.com/cpockrandt/genmap/wiki>`_.

For questions or feature requests feel free to open an issue on GitHub or send an e-mail to
``christopher.pockrandt [ÄT] fu-berlin.de``.

Christopher Pockrandt, Mai Alzamel, Costas S. Iliopoulos, Knut Reinert. **GenMap: Fast and Exact Computation of Genome Mappability**. `bioRxiv`_, presented on RECOMB-Seq, 2019.

.. _bioRxiv: https://doi.org/10.1101/611160

Binaries
^^^^^^^^

Your CPU must support the ``POPCNT`` instruction.
If you have a modern CPU, you can go with the optimized 64 bit version that additionally uses up to ``SSE4`` (MMX, SSE, SSE2, SSE3, SSSE3, SSE4).
This improves the running time by 10 %.
To verify whether your CPU supports these instructions sets you can check the output of
``cat /proc/cpuinfo | grep -E "mmx|sse|popcnt"`` (Linux) or
``sysctl -a | grep -i -E "mmx|sse|popcnt"`` (Mac).

.. Source of linux.svg: https://svgsilh.com/image/2025536.html
.. Source of apple.svg: https://svgsilh.com/image/2962084.html

+---------------------------------+---------------------------+------------------+-----------------------------+
| **Platform**                    | **Download**              | **Version**      | **Additional requirements** |
+---------------------------------+---------------------------+------------------+-----------------------------+
| .. image:: .github/linux.svg    | `Linux 64 bit`_           | 1.0 (2019-05-27) | \-                          |
+   :alt: Download Linux binaries +---------------------------+------------------+-----------------------------+
|   :height: 60px                 | `Linux 64 bit optimized`_ | 1.0 (2019-05-27) | requires up to SSE4         |
+---------------------------------+---------------------------+------------------+-----------------------------+
| .. image:: .github/apple.svg    | `Mac 64 bit`_             | 1.0 (2019-05-27) | \-                          |
+   :alt: Download Mac binaries   +---------------------------+------------------+-----------------------------+
|   :height: 60px                 | `Mac 64 bit optimized`_   | 1.0 (2019-05-27) | requires up to SSE4         |
+---------------------------------+---------------------------+------------------+-----------------------------+

.. _Linux 64 bit: https://github.com/cpockrandt/genmap/releases/download/genmap-v1.0/genmap-1.0-Linux-x86_64.zip
.. _Linux 64 bit optimized: https://github.com/cpockrandt/genmap/releases/download/genmap-v1.0/genmap-1.0-Linux-x86_64-sse4.zip
.. _Mac 64 bit: https://github.com/cpockrandt/genmap/releases/download/genmap-v1.0/genmap-1.0-Darwin-x86_64.zip
.. _Mac 64 bit optimized: https://github.com/cpockrandt/genmap/releases/download/genmap-v1.0/genmap-1.0-Darwin-x86_64-sse4.zip

Building from source
^^^^^^^^^^^^^^^^^^^^

Please note that building from source can easily take 10 minutes and longer depending on your machine and compiler.

::

    $ git clone --recursive https://github.com/cpockrandt/genmap.git
    $ mkdir genmap-build && cd genmap-build
    $ cmake ../genmap -DCMAKE_BUILD_TYPE=Release
    $ make genmap
    $ ./bin/genmap

If you are using a very old version of Git (< 1.6.5) the flag ``--recursive`` does not exist.
In this case you need to clone the submodule separately before you can run ``cmake``:

::

    $ git clone https://github.com/cpockrandt/genmap.git
    $ cd genmap
    $ git submodule update --init --recursive

Requirements
""""""""""""

Operating System
  GNU/Linux, Mac

Architecture
  Intel/AMD platforms that support ``POPCNT``

Compiler
  GCC ≥ 4.9, LLVM/Clang ≥ 3.8

Build system
  CMake ≥ 3.0

Language support
  C++14

Mappability example
^^^^^^^^^^^^^^^^^^^

Below you can see the (4,1)-mappability and frequency ``M`` and ``F`` of the nucleotide sequence ``T = ATCTAGCTTGCTAATCTA``.
Only mismatches (Hamming distance) are considered.
GenMap can also allow for insertions and deletions (Edit/Levenshtein distance, coming soon).

.. TODO: smaller example s.t. no scrolling is necessary

+----------+-------+-------+-------+-------+-------+-----+-----+-----+-----+-------+-------+-------+-------+------+-------+-------+-------+-------+
| **i**    |   0   |   1   |   2   |   3   |   4   |  5  |  6  |  7  |  8  |   9   |   10  |   11  |   12  |  13  |   14  |   15  |   16  |   17  |
+----------+-------+-------+-------+-------+-------+-----+-----+-----+-----+-------+-------+-------+-------+------+-------+-------+-------+-------+
| **T[i]** |   A   | **T** | **C** | **T** | **A** |  G  |  C  |  T  |  T  | **G** | **C** | **T** | **A** |   A  | **T** | **C** | **T** | **A** |
+----------+-------+-------+-------+-------+-------+-----+-----+-----+-----+-------+-------+-------+-------+------+-------+-------+-------+-------+
| **M[i]** |  0.33 |  0.33 |  0.33 |  0.5  |  0.25 | 0.5 | 0.5 | 0.5 | 0.5 |  0.25 |  0.5  |  1.0  |  1.0  | 0.33 |  0.33 |   0   |   0   |   0   |
+----------+-------+-------+-------+-------+-------+-----+-----+-----+-----+-------+-------+-------+-------+------+-------+-------+-------+-------+
| **F[i]** |   3   |   3   |   3   |   2   |   4   |  2  |  2  |  2  |  2  |   4   |   2   |   1   |   1   |   3  |   3   |   0   |   0   |   0   |
+----------+-------+-------+-------+-------+-------+-----+-----+-----+-----+-------+-------+-------+-------+------+-------+-------+-------+-------+

The mappability value ``M[1] = 0.33`` means that the 4-mer starting at position 1 ``T[1..3] = TCTA`` occurs three times in the sequence with up to one mismatch, namely at positions ``1 (TCTA)``, ``9 (GCTA)`` and ``14 (TCTA)``.

The mappability can be exported in various formats that allow post-processing or display in genome browsers.

Getting started
^^^^^^^^^^^^^^^

Building the index
""""""""""""""""""

At first you have to build an index of the fasta file(s) whose mappability you want to compute.
This step only has to performed once.
You might want to check out prebuilt indices for `download <#pre-built-indices>`_.

::

    $ ./genmap index -F /path/to/fasta.fasta -I /path/to/index/folder

A new folder ``/path/to/index/folder`` will be created to store the index and all associated files.

There are two algorithms that can be chosen for index construction.
One uses RAM (radix), one uses secondary memory (skew).
Depending on the quota and main memory limitations you can choose the appropriate algorithm with ``-A radix`` or
``-A skew``.
For skew you can change the location of the temp directory via the environment variable (e.g., to choose a directory
with more quota):

::

   $ export TMPDIR=/somewhere/else/with/more/space

Computing the mappability
"""""""""""""""""""""""""

To compute the (30,2)-mappability of the previously indexed genome, simply run:

::

    $ ./genmap map -E 2 -K 30 -I /path/to/index/folder -O /path/to/output/folder -t -w -b

This will create a ``text``, ``wig`` and ``bed`` file in ``/path/to/output/folder`` storing the computed mappability in
different formats. You can formats that are not required by omitting the corresponding flags ``-t`` ``-w`` or ``-b``.

Instead of the mappability, the frequency can be outputted, you only have to add the flag ``-fl`` to the previous
command.

Help pages and examples
"""""""""""""""""""""""

A detailed list of arguments and explanations can be retrieved with ``--help``:

::

    $ ./genmap --help
    $ ./genmap index --help
    $ ./genmap map --help

More detailed examples can be found in the Wiki.

Pre-built indices
^^^^^^^^^^^^^^^^^

Building an index on a large genome takes some time and requires a lot of space. Hence, we provide indexed genomes for download.
If you need other genomes indexed and do not have the computational resources, please send an e-mail to ``christopher.pockrandt [ÄT] fu-berlin.de``.

+------------------------------------------+-----------------------------+--------------------+
| **Genome**                               | **Index size (compressed)** | **Download**       |
+------------------------------------------+-----------------------------+--------------------+
| Human GRCh38 (`hg38 patch 13`_)          | 6.6 GB                      | `GRCh38 index`_    |
+------------------------------------------+-----------------------------+--------------------+
| Human GRCh37 (`hg19 patch 13`_)          | 6.4 GB                      | `GRCh37 index`_    |
+------------------------------------------+-----------------------------+--------------------+
| Mouse GRCm38 (`mm10 patch 6`_)           | 5.7 GB                      | `GRCm38 index`_    |
+------------------------------------------+-----------------------------+--------------------+
| Fruitfly D. melanogaster (`dm6 rel. 6`_) | 0.3 GB                      | `dm6 index`_       |
+------------------------------------------+-----------------------------+--------------------+
| Worm C. elegans (`ce11 WBcel235`_)       | 0.2 GB                      | `ce11 index`_      |
+------------------------------------------+-----------------------------+--------------------+

.. | Barley (`hordeum vulgare`_)              | x.x GB                      | `hv index`_        |
.. +------------------------------------------+-----------------------------+--------------------+

.. sequence: ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.28_GRCh38.p13/GCA_000001405.28_GRCh38.p13_genomic.fna.gz
.. _`hg38 patch 13`:   https://www.ncbi.nlm.nih.gov/assembly/GCA_000001405.28
.. sequence: ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_genomic.fna.gz
.. _`hg19 patch 13`:   https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.25
.. sequence: ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.26_GRCm38.p6/GCF_000001635.26_GRCm38.p6_genomic.fna.gz
.. _`mm10 patch 6`:    https://www.ncbi.nlm.nih.gov/assembly/GCF_000001635.26
.. sequence: ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/215/GCF_000001215.4_Release_6_plus_ISO1_MT/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna.gz
.. _`dm6 rel. 6`:      https://www.ncbi.nlm.nih.gov/assembly/GCF_000001215.4
.. sequence: ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/985/GCF_000002985.6_WBcel235/GCF_000002985.6_WBcel235_genomic.fna.gz
.. _`ce11 WBcel235`:   https://www.ncbi.nlm.nih.gov/assembly/GCF_000002985.6
.. sequence: ftp://ftp.ensemblgenomes.org/pub/plants/release-42/fasta/hordeum_vulgare/dna/Hordeum_vulgare.IBSC_v2.dna.toplevel.fa.gz
.. _`hordeum vulgare`: https://plants.ensembl.org/Hordeum_vulgare/Info/Index

.. _`GRCh38 index`: http://ftp.imp.fu-berlin.de/pub/cpockrandt/genmap/indices/grch38-dna5.tar.gz
.. _`GRCh37 index`: http://ftp.imp.fu-berlin.de/pub/cpockrandt/genmap/indices/grch37-dna5.tar.gz
.. _`GRCm38 index`: http://ftp.imp.fu-berlin.de/pub/cpockrandt/genmap/indices/mm10-dna5.tar.gz
.. _`dm6 index`:    http://ftp.imp.fu-berlin.de/pub/cpockrandt/genmap/indices/dm6-dna5.tar.gz
.. _`ce11 index`:   http://ftp.imp.fu-berlin.de/pub/cpockrandt/genmap/indices/ce11-dna5.tar.gz
.. _`hv index`:     http://ftp.imp.fu-berlin.de/pub/cpockrandt/genmap/indices/hv-dna5.tar.gz

Changelog
^^^^^^^^^

1.0 (2019-05-27)
""""""""""""""""

* faster computation of mappability for 3 and 4 errors
* included the reference to the paper (preprint)
* minor fixes in documentation and error messages

0.9 (2019-03-23)
""""""""""""""""

* preliminary version of GenMap released
