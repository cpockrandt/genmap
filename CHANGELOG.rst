Changelog
^^^^^^^^^

1.3.0 (2020-06-xx)
""""""""""""""""""

Features
--------

* **Faster and more efficient index creation (libdivsufsort)**, requires about 6n main memory (n is the size of the input fasta files) and no additional secondary memory (libdivsufsort is the default algorithm. Skew7 is still supported, radixsort was removed)

Fixes
-----

* indexing will search for fasta files recursively (subfolders were not considered before)
* updated paper information (replaced bioRxiv with OUP Bioinformatics)

1.2.0 (2020-02-21)
""""""""""""""""""

Features
--------

* **BREAKING CHANGE!** k-mers are searched on both strands by default. Can be turned off with ``--no-reverse-complement``
* bedgraph output (``*.bg``) replaces bed output (``*.bed``) (bed output is still supported with ``--bed`` but deprecated and removed from the help string)
* allow user to specify a filename with ``--output`` if only a single fasta file has been indexed (previously only the directory could be specified)

Fixes
-----

* truncate fasta identifiers after first space
* allow indexing fasta files with ``*.fas`` filename ending
* runtime speedup when BED file is provided for computation on a subset of the input (``--selection``)
* reduced progress output on terminal when processing multiple fasta files

1.1.0 (2019-11-17)
""""""""""""""""""

* compute mappability of selected regions using a bed file
* suppress 0 values in BED and WIG files
* change default algorithm for indexing to Skew

1.0.2 (2019-09-04)
""""""""""""""""""

* BED output format fixed (end position was off by one, i.e. closed interval instead of half-closed interval)

1.0.1 (2019-06-11)
""""""""""""""""""

* ``--frequency-small`` would output unreadable ascii characters in ``--txt`` format
* some typo fixes

1.0.0 (2019-05-27)
""""""""""""""""""

* faster computation of mappability for 3 and 4 errors
* included the reference to the paper (preprint)
* minor fixes in documentation and error messages

0.9.0 (2019-03-23)
""""""""""""""""""

* preliminary version of GenMap released
