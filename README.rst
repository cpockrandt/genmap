GenMap - Fast and Exact Computation of Genome Mappability |buildstatus|
-----------------------------------------------------------------------

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

The corresponding paper will be uploaded to biorxiv.org in mid-March.
Until then major design changes of the interface and minor changes to its specification are possible.

Binaries
^^^^^^^^

Your CPU must support the ``POPCNT`` instruction.
If you have a modern CPU, you can go with the optimized 64 bit version that additionally uses ``SSE4``.
To verify whether your CPU supports ``POPCNT`` and ``SSE4``, you can check the output of ``cat /proc/cpuinfo | grep 'popcnt\|sse4'``.

.. Source of download.svg: https://svgsilh.com/image/2203950.html

+---------------------------------+---------------------+----------------------------------+
| .. image:: .github/download.svg | `64 bit`_           | requires ``POPCNT``              |
+   :alt: Download binaries       +---------------------+----------------------------------+
|   :width: 56px                  | `64 bit optimized`_ | requires ``POPCNT`` and ``SSE4`` |
+---------------------------------+---------------------+----------------------------------+

.. _64 bit: http://ftp.imp.fu-berlin.de/pub/cpockrandt/genmap-0.9-Linux-x86_64.zip
.. _64 bit optimized: http://ftp.imp.fu-berlin.de/pub/cpockrandt/genmap-0.9-Linux-x86_64-sse4.zip

Building from source (currently for Linux only)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

NOTE: Building from source can take up to 10 minutes depending on your machine and compiler.

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

Below you can see the (4,1)-frequency ``F`` of the nucleotide sequence ``T = ATCTAGCTTGCTAATCTA``.
Only mismatches (Hamming distance) are considered.
GenMap can also allow for insertions and deletions (Edit distance).

+----------+---+-------+-------+-------+-------+---+---+---+---+-------+-------+-------+-------+----+-------+-------+-------+-------+
| **i**    | 0 |   1   |   2   |   3   |   4   | 5 | 6 | 7 | 8 |   9   |   10  |   11  |   12  | 13 |   14  |   15  |   16  |   17  |
+----------+---+-------+-------+-------+-------+---+---+---+---+-------+-------+-------+-------+----+-------+-------+-------+-------+
| **T[i]** | A | **T** | **C** | **T** | **A** | G | C | T | T | **G** | **C** | **T** | **A** |  A | **T** | **C** | **T** | **A** |
+----------+---+-------+-------+-------+-------+---+---+---+---+-------+-------+-------+-------+----+-------+-------+-------+-------+
| **F[i]** | 3 |   3   |   3   |   2   |   4   | 2 | 2 | 2 | 2 |   4   |   2   |   1   |   1   |  3 |   3   |   0   |   0   |   0   |
+----------+---+-------+-------+-------+-------+---+---+---+---+-------+-------+-------+-------+----+-------+-------+-------+-------+

The frequency value ``F[1] = 3`` means that the 4-mer starting at position 1 ``T[1..3] = TCTA`` occurs three times in the sequence with up to one mismatch, namely at positions ``1 (TCTA)``, ``9 (GCTA)`` and ``14 (TCTA)``.

Getting started
^^^^^^^^^^^^^^^

Building the index
""""""""""""""""""

At first you have to build an index of the fasta file(s) whose mappability you want to compute.
This step only has to performed once.

::

    $ ./genmap index -G /path/to/fasta.fasta -I /path/to/index/folder

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
different formats. You can remove not required formats by ommitting the corresponding flags ``-t`` ``-w`` or ``-b``.

Instead of the mappability, the frequency can be outputted, you only have to add the flag ``-fl`` to the previous
command.

Help pages and examples
"""""""""""""""""""""""

A detailed list of arguments and explanations can be retrieved using ``--help``:

::

    $ ./genmap --help
    $ ./genmap index --help
    $ ./genmap map --help

More detailed examples can be found in the Wiki.
