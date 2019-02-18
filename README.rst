GenMap - Fast and Exact Computation of Genome Mappability |buildstatus|
---------------------------------------------------------

.. |BUILDSTATUS| image:: https://travis-ci.org/cpockrandt/genmap.svg?branch=master
    :target: https://travis-ci.org/cpockrandt/genmap

GenMap is a tool to compute the mappability of nucleotide sequences.
In particular, it computes the (k,e)-frequency, i.e., how often each k-mer from the sequence occurs with up to e errors in the sequence itself.
The (k,e)-mappability is the inverse of the (k,e)-frequency.
Hence, a mappability value of 1 at position i indicates that the k-mer in the sequence at position i occurs only once in the sequence with up to e errors.
A low mappability value indicates that this k-mer belongs to a repetitive region.

A small example is listed below, for detailed examples such as marker sequence computation on multiple fasta files, please check out the GitHub wiki (coming soon).

Binaries
^^^^^^^^

Coming soon.

Building from source
^^^^^^^^^^^^^^^^^^^^

NOTE: Building from source can take up to 10 minutes depending on your machine.

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
  GCC ≥ 4.9, LLVM/Clang ≥ 3.9

Build system
  CMake ≥ 3.0

Language support
  C++14

Getting started
^^^^^^^^^^^^^^^

Building the index
""""""""""""""""""

At first you have to build an index of the fasta file(s) whose mappability you want to compute.
This step only has to performed once.

::

    $ ./genmap index -G /path/to/fasta.fasta -I /path/to/index/folder

There are two algorithms that can be chosen for index construction.
One uses RAM (radixsort), one uses secondary memory (skew).
Depending on the quota and main memory limitations you can choose the appropriate algorithm with ``-A radixsort`` or ``-A skew``.
For skew you can change the location of the temp directory via the environment variable (e.g., to choose a directory with more quota):

::

   $ export TMPDIR=/somewhere/else/with/more/space

Computing the mappability
"""""""""""""""""""""""""

::

    $ ./genmap map ...

Help pages and examples
"""""""""""""""""""""""

A detailed list of arguments and explanations can be retrieved using ``--help``:

::

    $ ./genmap --help
    $ ./genmap index --help
    $ ./genmap mappability --help

More detailed examples will be coming soon in the wiki.
