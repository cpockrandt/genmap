#include "../include/libdivsufsort/divsufsort.hpp"

namespace seqan
{
    template <typename sa_t>
    struct AlgoDivSufSortTag {};

    // total in brackets is for genomes <4GB and Dna5
    // 1. Packed text is in memory (seqan): (0.375n resp. 0.25n) (total: 0.375n)
    // 2. Copy text to c string: 1n (total: 1.375n)
    // 3. Compute SA with libdivsufsort: 4n resp. 8n (total: 5.375n)
    // 4. Delete c string (total: 4.375n)
    // 5. Compute CSA from SA: Xn + Yn bytes (total: ?), we should do this with External<> (TODO)
    // 6. Create BWT and bit vector indicating sentinels: 1.125n (total: 5.5n + CSA, since it has not been saved yet and is not External<> yet)
    // 7. Delete SA (total: 1.5n)
    // 8. Build auxiliary data structures for BWT / bit vector

  template <typename TAlphabet, typename TSeqNo, typename TSeqPos, typename sa_t, typename TConfig>
  inline bool indexCreate(Index<StringSet<String<TAlphabet, Packed<> >, Owner<ConcatDirect<SizeSpec_<TSeqNo, TSeqPos> > > >,
                                FMIndex<AlgoDivSufSortTag<sa_t>, TConfig> > & index,
                          FibreSALF)
    {
        typedef StringSet<String<TAlphabet, Packed<> >, Owner<ConcatDirect<SizeSpec_<TSeqNo, TSeqPos> > > > TText;
        typedef Index<TText, FMIndex<Nothing, TConfig> >                           TIndex;
        typedef typename Fibre<TIndex, FibreTempSA>::Type                          TSA;
        typedef typename Size<TSA>::Type                                           TSASize;

//        std::cout << "LIBDIVSUFSORT\n";

//	      time_t tt;

        TText const & text = indexText(index);

        if (empty(text))
            return false;

        // compute number of sequences, cumulative sequence lengths, size of CSA, etc.
        uint64_t const nbr_sequences = length(text);
        std::vector<uint64_t> cum_seq_lengths(nbr_sequences + 1);
        uint64_t seq_id = 0;
        uint64_t csa_size = 0;
        for (auto const len : stringSetLimits(text))
        {
            cum_seq_lengths[seq_id] = len + seq_id; // stringSetLimits are cumulative values but do not count the sentinels
            if (seq_id > 0) // first entry is 0
            {
                uint64_t const textlength = cum_seq_lengths[seq_id] - cum_seq_lengths[seq_id - 1] - 1; // exclude sentinel
                csa_size += ((textlength - 1) / TConfig::SAMPLING) + 1; // == ceil(textlength / TConfig<>::SAMPLING)
            }
            ++seq_id;
        }
        sa_t const sequences_length_with_sentinels = cum_seq_lengths.back();

        // copy text to c string (with sentinels)
//        tt = time(NULL); printf("\n%s\tCreate C string text", ctime(&tt));
        uint8_t * ctext = static_cast<uint8_t *>(malloc(sizeof(uint8_t) * sequences_length_with_sentinels));

        for (uint64_t i = 0, j = 0; j < nbr_sequences; ++j)
        {
            uint64_t const seq_length = cum_seq_lengths[j + 1] - cum_seq_lengths[j] - 1; // -1 because we don't count the sentinel
            for (uint64_t k = 0; k < seq_length; ++k, ++i)
            {
                ctext[i] = ordValue(text[j][k]) + 1;
            }
            ctext[i] = 0; // sentinel
            ++i;
        }

        // compute full suffix array with libdivsufsort
//        tt = time(NULL); printf("\n%s\tBuild full SA", ctime(&tt));
        sa_t * sa = static_cast<sa_t *>(malloc(sizeof(sa_t) * sequences_length_with_sentinels));
        sdsl::divsufsort(ctext, sa, sequences_length_with_sentinels);
        // clear c string of text
        ::free(ctext);

        // Set the FMIndex LF as the CompressedSA LF.
        setFibre(indexSA(index), indexLF(index), FibreLF());

        // Create the compressed SA.
        // former: createCompressedSa(indexSA(index), tempSA, nbr_sequences);
//        tt = time(NULL); printf("\n%s\tBuild CSA", ctime(&tt));
        {
            typedef CompressedSA<TText, Nothing, TConfig>                  TCompressedSA;
            typedef typename Fibre<TCompressedSA, FibreSparseString>::Type TSparseSA;
            typedef typename Fibre<TSparseSA, FibreIndicators>::Type       TIndicators;
            typedef typename Fibre<TSparseSA, FibreValues>::Type           TValues;

            auto & compressedSA = indexSA(index);

            TSparseSA & sparseString = getFibre(compressedSA, FibreSparseString());
            TIndicators & indicators = getFibre(sparseString, FibreIndicators());
            TValues & values = getFibre(sparseString, FibreValues());

            resize(compressedSA, sequences_length_with_sentinels, Exact()); // resizes only indicators, not values.

            TSASize pos = 0;
            for (; pos < nbr_sequences; ++pos)
                setValue(indicators, pos, false);

            resize(values, csa_size);

            for (uint64_t counter = 0; pos < static_cast<uint64_t>(sequences_length_with_sentinels); ++pos)
            {
                auto const u = upper_bound(cum_seq_lengths.begin(), cum_seq_lengths.end(), sa[pos]);
                if (static_cast<uint64_t>(sa[pos]) + 1 != *u) // ignore sentinel positions
                {
                    TSeqNo const i1 = difference(cum_seq_lengths.begin(), u-1);
                    TSeqPos const i2 = sa[pos] - *(u-1);

                    if (i2 % TConfig::SAMPLING == 0)
                    {
                        assignValue(values, counter, Pair<TSeqNo, TSeqPos>(i1, i2));
                        setValue(indicators, pos, true);
                        ++counter;
                    }
                    else
                        setValue(indicators, pos, false);
                }
            }

//            tt = time(NULL); printf("\n%s\tUpdate CSA ranks", ctime(&tt));
            updateRanks(indicators);

            // TODO: test this with different sampling rates
            if (getRank(indicators, length(sparseString) - 1) != length(values))
            {
                std::cerr << "ERROR: It seems that the size of `values` has been precomputed incorrectly!\n";
                exit(12);
            }
        }

        // Create the LF table.
        // former: createLF(indexLF(index), text, tempSA);
//        tt = time(NULL); printf("\n%s\tBuild BWT", ctime(&tt));
        {
            auto & lf = indexLF(index);

            typedef LF<TText, Nothing, TConfig> TLF;
            typedef typename Value<TLF>::Type   TValue;
            typedef typename Size<TLF>::Type    TSize;

            // Clear assuming undefined state.
            clear(lf);

            // Compute prefix sum.
            prefixSums<TValue>(lf.sums, text);

            // Choose the sentinel substitute.
            _setSentinelSubstitute(lf);

            // Create and index BWT bwt for rank queries.
            // former: createRankDictionary(lf, text, sa);
            {
                // Resize the RankDictionary.
                resize(lf.sentinels, sequences_length_with_sentinels, Exact());
                resize(lf.bwt, sequences_length_with_sentinels, Exact()); // TODO: make sure that this does not allocate memory for precomputed ranks yet

                // Fill the sentinel positions (they are all at the beginning of the bwt).
                TSize i = 0;
                for (; i < nbr_sequences; ++i)
                {
                    // if (length(text[nbr_sequences - (i + 1)]) > 0) // not necessary, genmap removes empty sequences beforehand
                    // {
                        auto const u = upper_bound(cum_seq_lengths.begin(), cum_seq_lengths.end(), sa[i]);
                        uint64_t const i1 = difference(cum_seq_lengths.begin(), u-1);

                        setValue(lf.bwt, i, back(text[i1]));
                        setValue(lf.sentinels, i, false);
                    // }
                }

                // Compute the rest of the BWT
                for (; i < static_cast<uint64_t>(sequences_length_with_sentinels); ++i)
                {
                    auto const u = upper_bound(cum_seq_lengths.begin(), cum_seq_lengths.end(), sa[i]);
                    uint64_t const i1 = difference(cum_seq_lengths.begin(), u-1);
                    uint64_t const i2 = sa[i] - *(u-1);

                    if (i2 != 0)
                    {
                        setValue(lf.bwt, i, text[i1][i2 - 1]);
                        setValue(lf.sentinels, i, false);
                    }
                    else
                    {
                        setValue(lf.bwt, i, lf.sentinelSubstitute);
                        setValue(lf.sentinels, i, true);
                    }
                }
//                tt = time(NULL); printf("\n%s\tUpdate ranks for BWT", ctime(&tt));

                // Delete full suffix array
                ::free(sa);

                // Update all ranks.
                updateRanks(lf.bwt);
                // Update the auxiliary RankDictionary of sentinel positions.
                updateRanks(lf.sentinels);
            }

            // Add sentinels to prefix sum.
            for (TSize i = 0; i < length(lf.sums); ++i)
                lf.sums[i] += nbr_sequences;
        }

        return true;
    }

}
