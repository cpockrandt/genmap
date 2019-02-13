// ==========================================================================
//                                  lambda
// ==========================================================================
// Copyright (c) 2013-2017, Hannes Hauswedell <h2 @ fsfe.org>
// Copyright (c) 2016-2017, Knut Reinert and Freie Universität Berlin
// All rights reserved.
//
// This file is part of Lambda.
//
// Lambda is Free Software: you can redistribute it and/or modify it
// under the terms found in the LICENSE[.md|.rst] file distributed
// together with this file.
//
// Lambda is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//
// ==========================================================================
// lambda_indexer_misc.hpp: misc stuff for indexer
// ==========================================================================

#ifndef LAMBDA_INDEXER_MISC_HPP_
#define LAMBDA_INDEXER_MISC_HPP_

using namespace seqan;

// ============================================================================
// Parallel BWT construction
// ============================================================================

template <typename TText, typename TSSetSpec, typename TSpec, typename TConfig,
          typename TOtherText,
          typename TSA,
          typename TCallback>
void
createRankDictionaryProgress(LF<StringSet<TText, TSSetSpec>, TSpec, TConfig> & lf,
                             TOtherText const & text,
                             TSA const & sa,
                             TCallback && progress)
{
    typedef typename Value<TSA>::Type                       TSAValue;
    typedef typename Size<TSA>::Type                        TSize;

    // Resize the RankDictionary.
    TSize seqNum = countSequences(text);
    TSize totalLen = lengthSum(text);
    resize(lf.sentinels, seqNum + totalLen, Exact());
    resize(lf.bwt, seqNum + totalLen, Exact());

    // Fill the sentinel positions (they are all at the beginning of the bwt).
    for (TSize i = 0; i < seqNum; ++i)
    {
        if (length(text[seqNum - (i + 1)]) > 0)
        {
            setValue(lf.bwt, i, back(text[seqNum - (i + 1)]));
            setValue(lf.sentinels, i, false);
        }
    }

    /* Compute the rest of the bwt.*/

    // align the chunk_size to underlying word boundaries to prevent parallel write to word spanning chunk boundary
    uint64_t const chunkSize = std::max(static_cast<uint64_t>((length(sa) / omp_get_max_threads() / 64ull) * 64ull), uint64_t{1});
    uint64_t const twoPercent = std::max(chunkSize / 50, uint64_t{1});
    // the 0th thread might get an additional chunk because of the above alignment so we count from the 1st instead
    uint32_t const countThreadID = omp_get_max_threads() > 1 ? 1 : 0;

    SEQAN_OMP_PRAGMA(parallel for schedule(static, chunkSize))
    for (TSize i = 0; i < length(sa); ++i)
    {
        TSAValue pos;    // = SA[i];
        posLocalize(pos, sa[i], stringSetLimits(text));

        if (getSeqOffset(pos) != 0)
        {
            setValue(lf.bwt, i + seqNum, getValue(getValue(text, getSeqNo(pos)), getSeqOffset(pos) - 1));
            setValue(lf.sentinels, i + seqNum, false);
        }
        else
        {
            setValue(lf.bwt, i + seqNum, lf.sentinelSubstitute);
            setValue(lf.sentinels, i + seqNum, true);
        }

        if (((static_cast<uint32_t>(omp_get_thread_num()) == countThreadID) && ((i % chunkSize) % twoPercent == 0)))
            progress(((i % chunkSize) / twoPercent) * 2);
    }

   // Update all ranks.
   updateRanks(lf.bwt);
   // Update the auxiliary RankDictionary of sentinel positions.
   updateRanks(lf.sentinels);
}

template <typename TText, typename TSpec, typename TConfig, typename TOtherText, typename TSA, typename TCallback>
void
createLFProgress(LF<TText, TSpec, TConfig> & lf, TOtherText const & text, TSA const & sa, TCallback && progress)
{
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
    createRankDictionaryProgress(lf, text, sa, progress);

    // Add sentinels to prefix sum.
    TSize sentinelsCount = countSequences(text);
    for (TSize i = 0; i < length(lf.sums); ++i)
        lf.sums[i] += sentinelsCount;

    progress(100);
}

#endif // LAMBDA_INDEXER_MISC_HPP_
