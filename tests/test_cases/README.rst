GenMap - Fast and Exact Computation of Genome Mappability
---------------------------------------------------------

Details on test cases
^^^^^^^^^^^^^^^^^^^^^

There are a few hand-written test cases to check the output formats of GenMap and possible edge cases.

Single fasta file with a single sequence (i.e., only one chromosome)
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

case 1a
  3-mers with 0 errors without the reverse complement

case 1b
  3-mers with 0 errors with the reverse complement

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

Single fasta file in directory with multiple sequences
""""""""""""""""""""""""""""""""""""""""""""""""""""""

Coming soon.

Multiple fasta files in directory with single and multiple sequences
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Coming soon (including and exluding ``--exclude-pseudo``).

.. [1] This is important due to the design of the algorithm. All sequences in a fasta file are concatenated and each k-mer in the concatenated string is searched. False positives are later filtered out, i.e., windows that spanned multiple sequences.
