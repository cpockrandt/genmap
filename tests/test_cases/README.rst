GenMap - Fast and Exact Computation of Genome Mappability
---------------------------------------------------------

Details on test cases
^^^^^^^^^^^^^^^^^^^^^

There are a few hand-written test cases to check the output formats of GenMap and possible edge cases. If not stated otherwise, the alphabet is Dna4 (A, C, G, T).

Single fasta file with a single sequence (i.e., only one chromosome)
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

case 1a
  3-mers with 0 errors without the reverse complement

case 1b
  3-mers with 0 errors with the reverse complement

case 1c
  3-mers with 0 errors without the reverse complement. **Dna5** alphabet (including Ns).

case 1d
  3-mers with 0 errors with the reverse complement. **Dna5** alphabet (including Ns).

case 1e
  3-mers with 1 error without the reverse complement. **Dna5** alphabet (including Ns).

case 1f
  3-mers with 1 error with the reverse complement. **Dna5** alphabet (including Ns).

case 1g
  3-mers with 1 error with the reverse complement. **Dna5** alphabet (including Ns). Computing a subset defined in a bed file.

Single fasta file with multiple sequences
"""""""""""""""""""""""""""""""""""""""""

case 2a
  4-mers with 0 errors without the reverse complement. k-mers spanning two sequences occur within the sequence itself [1]_.

case 2b
  case 2a with reverse complement

case 2c
  4-mers with 0 errors without the reverse complement. Contains sequences that are shorter or equal to k. The last sequence is also shorter than k. k-mers spanning two sequences occur within some sequence, as well as in some sequences reverse complement [1]_.

case 2d
  case 2c with reverse complement

case 2e
  case 2d but only computing a subset defined in a bed file.

Multiple fasta files in directory with single and multiple sequences
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

case 3a
    4-mers with 0 errors without reverse complement. Multiple fasta files. Some only have one sequence (chromosome), some have multiple. There are k-mers that
      - occur once only in a single sequence in only one fasta file (TTTT)
      - occur multiple times only in a single sequence in only one fasta file (ACGT)
      - occur multiple times in different sequences of only one fasta file (CGTA)
      - occur multiple times in different sequences of multiple fasta files (at most one hit per fasta file) (AGGA)
      - occur multiple times in different sequences of multiple fasta files (multiple hits per fasta file) (ACCA, AAGG)

case 3b
    case 3a with reverse complement

case 3c
    case 3a with ``--exclude-pseudo``

case 3d
    case 3c with reverse complement

case 3e
    case 3d but only computing a subset defined in a bed file.

Single fasta file in directory with multiple sequences
""""""""""""""""""""""""""""""""""""""""""""""""""""""

Coming soon (including and exluding ``--exclude-pseudo``).

.. [1] This is important due to the design of the algorithm. All sequences in a fasta file are concatenated and each k-mer in the concatenated string is searched. False positives are later filtered out, i.e., windows that spanned multiple sequences.
