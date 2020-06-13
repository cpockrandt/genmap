#include "../include/divsufsort.hpp"

namespace seqan
{
    struct AlgoDivSufSortTag {};

    // total in brackets is for genomes <4GB and Dna5 and sampling rate of 10
    // 1. Load text into memory (seqan, unpacked): 1n or less (0.5n resp. 0.25n) (total: 1n)
    // 2. Copy text to c string: 1n (total: 2n)
    // 2b. Delete text from memory (seqan): (total: 1n)
    // 3. Compute SA with libdivsufsort: 4n or 8n (total: 5n)
    // 4. Compute CSA from SA: Xn + Yn bytes, possible as external string(?) (total: 5n)
    // 4. compute BWT (with c string and SA): 1n (total: 6n)
    // 5a. delete c text and SA
    // 5b. build dictionary on BWT

    // try to use text instead of ctext (if we pack it, we can free ctext instead!
    // also we don't have to store/load the text for bidirectional indices
    // but then we will have 1n + 1n + (4n or 8n) peak instead of 1n less!
    // let's free ctext after SA construction and see whether packing text gives us a significant performance drop

    template <typename TText, typename TLengthSum, unsigned LEVELS, unsigned WORDS_PER_BLOCK>
    inline bool indexCreate(Index<TText, FMIndex<AlgoDivSufSortTag, GemMapFastFMIndexConfig<void, TLengthSum, LEVELS, WORDS_PER_BLOCK> > > & index, FibreSALF)
    {
        typedef Nothing TSpec;
        typedef GemMapFastFMIndexConfig<void, TLengthSum, LEVELS, WORDS_PER_BLOCK> TConfig;

        typedef std::conditional_t<std::is_same<TLengthSum, uint32_t>::value, int32_t, int64_t> sa_t;

	time_t tt;

        std::cout << "Using libdivsufsort\n";

        typedef Index<TText, FMIndex<TSpec, TConfig> >  TIndex;
        typedef typename Fibre<TIndex, FibreTempSA>::Type            TTempSA;

        TText const & text = indexText(index);

        if (empty(text))
            return false;

        // compute number of sequences, cumulative sequence lengths, etc.
        uint64_t const nbr_sequences = length(text);
        std::vector<uint64_t> cum_seq_lengths(nbr_sequences + 1);
        uint64_t seq_id = 0;
        uint64_t csa_size = 0;
        for (auto & len : stringSetLimits(text))
        {
            cum_seq_lengths[seq_id] = len + seq_id; // stringSetLimits are cumulative values but do not count the sentinels
//            std::cout << cum_seq_lengths[seq_id] << '\n';
            if (seq_id > 0) // first entry is 0
            {
                uint64_t const textlength = cum_seq_lengths[seq_id] - cum_seq_lengths[seq_id - 1] - 1; // exclude sentinel
                csa_size += ((textlength - 1) / TConfig::SAMPLING) + 1; // == ceil(textlength / TConfig<>::SAMPLING)
            }
            ++seq_id;
        }
        sa_t const sequences_length_with_sentinels = cum_seq_lengths.back();
//        uint64_t const sequences_length_without_sentinels = cum_seq_lengths.back() - nbr_sequences; // length of all sequences without sentinels

//        std::cout << sequences_length_without_sentinels << ' ' << lengthSum(text) << '\n';

//        int64_t len = 0;
//        lengths.reserve(length(text) + 1);
//        lengths.push_back(0);
//    std::cout << "Lengths: 0";
//        for (uint64_t j = 0; j < length(text); ++j)
//        {
//            len += length(text[j]) + 1; // for sentinel
//            lengths.push_back(len);
//        std::cout << ' ' << len;
//        }

//        std::cout << '\n';
//
//        std::cout << "Limits: ";
//        for (auto & x : stringSetLimits(text))
//            std::cout << x << ' ';
//        std::cout << '\n';
//
//        exit (17);

        // get text in c string (with sentinels)
        tt = time(NULL); printf("\n%s\tCreate C string text\n", ctime(&tt));
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
        tt = time(NULL); printf("\n%s\tBuild full SA\n", ctime(&tt));
        sa_t * sa = static_cast<sa_t *>(malloc(sizeof(sa_t) * sequences_length_with_sentinels));
        sdsl::divsufsort(ctext, sa, sequences_length_with_sentinels);
        // TODO: clear c string of text
        ::free(ctext);

//    std::cout << "SA(lib):";
//    for (int64_t i = 0; i < sequences_length_with_sentinels; ++i)
//        std::cout << ' ' << sa[i];
//    std::cout << '\n';

        // Set the FMIndex LF as the CompressedSA LF.
        setFibre(indexSA(index), indexLF(index), FibreLF());

        // copy SA into tempSA (not needed later, just for later debugging)
//        TTempSA tempSA;
//        resize(tempSA, lengthSum(text), Exact());
//        //createSuffixArray(tempSA, text, TAlgo());
//        uint64_t pos = 0;

//std::cout << "SA(new): ";
//        for (int64_t i = 0; i < sequences_length_with_sentinels; ++i)
//        {
//            auto const u = upper_bound(cum_seq_lengths.begin(), cum_seq_lengths.end(), sa[i]);
//            uint64_t const i1 = difference(cum_seq_lengths.begin(), u-1);
//            uint64_t const i2 = sa[i] - *(u-1);
//            if (static_cast<uint64_t>(sa[i]) + 1 != *u)
//            {
////                tempSA[pos] = Pair<typename Size<TIndex>::Type>({i1, i2});
//            std::cout << " (" << i1 << ", " << i2 << ")";
//            }
//            else
//            {
//                std::cout << " (" << i1 << ", " << i2 << ")*";
//            }
//        }
//        std::cout << '\n';

        // Create the compressed SA.
        // former: createCompressedSa(indexSA(index), tempSA, nbr_sequences);
        tt = time(NULL); printf("\n%s\tBuild CSA\n", ctime(&tt));
        {
            typedef CompressedSA<TText, TSpec, TConfig>        TCompressedSA;
            typedef typename Size<TTempSA>::Type                                TSASize;
            typedef typename Fibre<TCompressedSA, FibreSparseString>::Type  TSparseSA;
            typedef typename Fibre<TSparseSA, FibreIndicators>::Type        TIndicators;
            typedef typename Fibre<TSparseSA, FibreValues>::Type            TValues;

            auto & compressedSA = indexSA(index);

            TSparseSA & sparseString = getFibre(compressedSA, FibreSparseString());
            TIndicators & indicators = getFibre(sparseString, FibreIndicators());
            TValues & values = getFibre(sparseString, FibreValues());

            resize(compressedSA, sequences_length_with_sentinels, Exact()); // resizes only indicators, not values. textlength includes sentinels

            TSASize pos = 0;
            for (; pos < nbr_sequences; ++pos)
                setValue(indicators, pos, false);

            resize(values, csa_size);

            for (uint64_t counter = 0; pos < static_cast<uint64_t>(sequences_length_with_sentinels); ++pos)
            {
                auto const u = upper_bound(cum_seq_lengths.begin(), cum_seq_lengths.end(), sa[pos]);
                if (static_cast<uint64_t>(sa[pos]) + 1 != *u) // ignore sentinel positions
                {
                    uint64_t const i1 = difference(cum_seq_lengths.begin(), u-1);
                    uint64_t const i2 = sa[pos] - *(u-1);

                    if (i2 % TConfig::SAMPLING == 0)
                    {
                        typedef typename GetValue<TValues>::Type MyPair;
                        MyPair p({i1, i2});
                        assignValue(values, counter, p);
                        setValue(indicators, pos, true);
                        ++counter;
                    }
                    else
                        setValue(indicators, pos, false);
                }
            }

            tt = time(NULL); printf("\n%s\tUpdate CSA ranks\n", ctime(&tt));
            updateRanks(indicators);

            if (getRank(indicators, length(sparseString) - 1) != length(values))
            {
                std::cerr << "ERROR: It seems that the size of `values` has been precomputed incorrectly!\n";
                exit(12);
            }
        }

        //std::cout << "SA : ";
        //for (uint64_t i = 0; i < length(tempSA); ++i)
        //    std::cout << "(" << tempSA[i].i1 << ", " << tempSA[i].i2 << ") ";
        //std::cout << '\n';

        // Create the LF table.
        // former: createLF(indexLF(index), text, tempSA);
        tt = time(NULL); printf("\n%s\tBuild BWT\n", ctime(&tt));
        {
            auto & lf = indexLF(index);

            typedef LF<TText, TSpec, TConfig>                          TLF;
            typedef typename Value<TLF>::Type                          TValue;
            typedef typename Size<TLF>::Type                           TSize;

            // Clear assuming undefined state.
            clear(lf);

            // Compute prefix sum.
            prefixSums<TValue>(lf.sums, text);

            // Choose the sentinel substitute.
            _setSentinelSubstitute(lf);

            // Create and index BWT bwt for rank queries.
            // former: createRankDictionary(lf, text, sa);
            {
                typedef typename Size<TTempSA>::Type                        TSize;

                // Resize the RankDictionary.
                resize(lf.sentinels, sequences_length_with_sentinels, Exact());
                resize(lf.bwt, sequences_length_with_sentinels, Exact());

                // Fill the sentinel positions (they are all at the beginning of the bwt).
                TSize i = 0;
                for (; i < nbr_sequences; ++i)
                {
                    // if (length(text[nbr_sequences - (i + 1)]) > 0) // TODO: not necessary
                    // {
                        auto const u = upper_bound(cum_seq_lengths.begin(), cum_seq_lengths.end(), sa[i]);
                        uint64_t const i1 = difference(cum_seq_lengths.begin(), u-1);

                        setValue(lf.bwt, i, back(text[i1]));
                        setValue(lf.sentinels, i, false);
                    // }
                }

                // Compute the rest of the bwt.
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
                tt = time(NULL); printf("\n%s\tUpdate ranks for WT\n", ctime(&tt));
                // Update all ranks.
                updateRanks(lf.bwt);
                // Update the auxiliary RankDictionary of sentinel positions.
                updateRanks(lf.sentinels);
            }

//            std::cout << "BWT: ";
//            uint64_t i = 0;
//            for (; i < sequences_length_with_sentinels; ++i)
//                std::cout << getValue(lf.bwt, i);
//            std::cout << '\n';

            // Add sentinels to prefix sum.
            for (TSize i = 0; i < length(lf.sums); ++i)
                lf.sums[i] += nbr_sequences;
        }



        ::free(sa);


        return true;
    }

}
