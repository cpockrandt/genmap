GenMap |VERSION|: Fast and Exact Computation of Genome Mappability
==================================================================

.. image:: https://img.shields.io/conda/dn/bioconda/genmap.svg?style=flag&label=BioConda%20install
    :target: https://anaconda.org/bioconda/genmap
    :alt: BioConda Install
.. image:: https://img.shields.io/github/downloads/cpockrandt/genmap/total.svg
    :target: https://github.com/cpockrandt/genmap/releases/latest
    :alt: Github All Releases
.. image:: https://travis-ci.org/cpockrandt/genmap.svg?branch=master
    :target: https://travis-ci.org/cpockrandt/genmap
    :alt: Travis CI
.. image:: https://img.shields.io/badge/License-BSD%203--Clause-blue.svg
    :target: https://opensource.org/licenses/BSD-3-Clause
    :alt: BSD3 License

.. contents::
   :local:
   :depth: 2

Introduction
^^^^^^^^^^^^

GenMap computes the uniqueness of k-mers for each position in the genome while allowing for up to e mismatches.
More formally, the uniqueness or (k,e)-mappability can be described for every position as the reciprocal value of how often each k-mer occurs approximately in the genome, i.e., with up to e mismatches.
Hence, a mappability value of 1 at position i indicates that the k-mer in the sequence at position i occurs only once in the sequence with up to e errors.
A low mappability value indicates that this k-mer belongs to a repetitive region.
GenMap can be applied to single or multiple genomes and helps finding regions that are unique or shared by many or all genomes.

Below you can see the (4,1)-mappability and frequency ``M`` and ``F`` of the nucleotide sequence ``T = ATCTAGGCTAATCTA``.
The mappability value ``M[1] = 0.33`` means that the 4-mer starting at position 1 ``T[1..3] = TCTA`` occurs three times in the sequence with up to one mismatch: at positions ``1 (TCTA)``, ``6 (GCTA)`` and ``11 (TCTA)``.

.. image:: .github/example.png
   :align: center
   :alt: example of mappability

The mappability can be exported in various formats that allow post-processing or display in genome browsers.
A small example on how to run GenMap is listed below, further details are on the `GitHub Wiki pages <https://github.com/cpockrandt/genmap/wiki>`_.
For questions or feature requests feel free to open an issue on GitHub or send an e-mail to ``christopher.pockrandt [ÄT] fu-berlin.de``.

Christopher Pockrandt, Mai Alzamel, Costas S. Iliopoulos, Knut Reinert. **GenMap: Fast and Exact Computation of Genome Mappability**. `bioRxiv`_, presented on RECOMB-Seq, 2019.

.. _bioRxiv: https://doi.org/10.1101/611160

Installation
^^^^^^^^^^^^

Bioconda
""""""""

::

    $ conda install -c bioconda genmap

Binaries
""""""""

Your CPU must support the ``POPCNT`` instruction.
If you have a modern CPU, you can go with the optimized 64 bit version that additionally uses SSE4.
This improves the running time by 10 %.
To verify whether your CPU supports these instructions sets you can check the output of
``cat /proc/cpuinfo | grep -E "mmx|sse|popcnt"`` (Linux) or
``sysctl -a | grep -i -E "mmx|sse|popcnt"`` (Mac).

.. Source of linux.svg: https://svgsilh.com/image/2025536.html
.. Source of apple.svg: https://svgsilh.com/image/2962084.html

+---------------------------------+---------------------------+--------------------------+-----------------------------+
| **Platform**                    | **Download**              | **Version**              | **Additional requirements** |
+---------------------------------+---------------------------+--------------------------+-----------------------------+
| .. image:: .github/linux.svg    | `Linux 64 bit`_           | |VERSION| (|BUILD_DATE|) | \-                          |
+   :alt: Download Linux binaries +---------------------------+--------------------------+-----------------------------+
|   :height: 60px                 | `Linux 64 bit optimized`_ | |VERSION| (|BUILD_DATE|) | requires SSE4               |
+---------------------------------+---------------------------+--------------------------+-----------------------------+
| .. image:: .github/apple.svg    | `Mac 64 bit`_             | |VERSION| (|BUILD_DATE|) | \-                          |
+   :alt: Download Mac binaries   +---------------------------+--------------------------+-----------------------------+
|   :height: 60px                 | `Mac 64 bit optimized`_   | |VERSION| (|BUILD_DATE|) | requires SSE4               |
+---------------------------------+---------------------------+--------------------------+-----------------------------+

.. _Linux 64 bit: https://github.com/cpockrandt/genmap/releases/download/genmap-v1.1.0/genmap-1.1.0-Linux-x86_64.zip
.. _Linux 64 bit optimized: https://github.com/cpockrandt/genmap/releases/download/genmap-v1.1.0/genmap-1.1.0-Linux-x86_64-sse4.zip
.. _Mac 64 bit: https://github.com/cpockrandt/genmap/releases/download/genmap-v1.1.0/genmap-1.1.0-Darwin-x86_64.zip
.. _Mac 64 bit optimized: https://github.com/cpockrandt/genmap/releases/download/genmap-v1.1.0/genmap-1.1.0-Darwin-x86_64-sse4.zip

.. |VERSION| replace:: 1.1.0
.. |BUILD_DATE| replace:: 2019-11-17

Building from source
""""""""""""""""""""

Please note that building from source can easily take 10 minutes and longer depending on your machine and compiler.

::

    $ git clone --recursive https://github.com/cpockrandt/genmap.git
    $ mkdir genmap-build && cd genmap-build
    $ cmake ../genmap -DCMAKE_BUILD_TYPE=Release
    $ make genmap

You can install genmap as follows

::

    $ sudo make install
    $ genmap

or run the binary directly:

::

    $ ./genmap

If you are using a very old version of Git (< 1.6.5) the flag ``--recursive`` does not exist.
In this case you need to clone the submodule separately before you can run ``cmake``:

::

    $ git clone https://github.com/cpockrandt/genmap.git
    $ cd genmap
    $ git submodule update --init --recursive

**Requirements**

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

Getting started
^^^^^^^^^^^^^^^

Building the index
""""""""""""""""""

At first you have to build an index of the fasta file(s) whose mappability you want to compute.
This step only has to be performed once.
You might want to check out pre-built indices available for `download <#pre-built-indices>`_.

::

    $ ./genmap index -F /path/to/fasta.fasta -I /path/to/index/folder

A new folder ``/path/to/index/folder`` will be created to store the index and all associated files.

There are two algorithms that can be chosen for index construction.
One uses RAM (radix), one uses secondary memory (skew).
Depending on the quota and main memory limitations you can choose the appropriate algorithm with ``-A radix`` or
``-A skew``.
It is recommended to use Skew, since Radix is comparison-based and therefore significantly slower on repetitive data.
For skew you can change the location of the temp directory via the environment variable (e.g., to choose a directory
with more quota):

::

   $ export TMPDIR=/somewhere/else/with/more/space

Computing the mappability
"""""""""""""""""""""""""

To compute the (30,2)-mappability of the previously indexed genome, simply run:

::

    $ ./genmap map -K 30 -E 2 -I /path/to/index/folder -O /path/to/output/folder -t -w -b

This will create a ``text``, ``wig`` and ``bed`` file in ``/path/to/output/folder`` storing the computed mappability in
different formats. You can omit formats that are not required by removing the corresponding flags ``-t`` ``-w`` or ``-b``.

Instead of the mappability, the frequency can be outputted, you only have to add the flag ``-fl`` to the previous
command.

Help pages and examples
"""""""""""""""""""""""

A detailed list of arguments and explanations can be retrieved with ``--help``:

::

    $ ./genmap --help
    $ ./genmap index --help
    $ ./genmap map --help

More detailed examples can be found in the `Wiki <https://github.com/cpockrandt/genmap/wiki>`_.

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
