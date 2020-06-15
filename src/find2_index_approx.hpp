// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2018, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of NVIDIA Corporation nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL NVIDIA CORPORATION BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Christopher Pockrandt <github@cpockrandt.de>
// ==========================================================================
// Modified Optimum Search Schemes for GenMap (GM)
// ==========================================================================

#pragma once

namespace seqan {

struct OptimalSearchDynGM
{
    std::vector<uint8_t> pi; // order of the blocks. permutation of [1..n]
    std::vector<uint8_t> l; // minimum number of errors at the end of the corresponding block
    std::vector<uint8_t> u; // maximum number of errors at the end of the corresponding block

    std::vector<uint32_t> blocklength; // cumulated values / prefix sums
    uint32_t startPos;
};

// std::vector<OptimalSearchDynGM<1> > VALUE { {{ {{1}}, {{0}}, {{10}}, {{0}}, 0 }} };

template <size_t N>
struct OptimalSearchGM
{
    std::array<uint8_t, N> pi; // order of the blocks. permutation of [1..n]
    std::array<uint8_t, N> l; // minimum number of errors at the end of the corresponding block
    std::array<uint8_t, N> u; // maximum number of errors at the end of the corresponding block

    std::array<uint32_t, N> blocklength; // cumulated values / prefix sums
    uint32_t startPos;
};

template <size_t max, typename TVoidType = void>
struct OptimalSearchSchemesGM;

template <typename TVoidType>
struct OptimalSearchSchemesGM<0, TVoidType>
{
    static constexpr std::array<OptimalSearchGM<1>, 1> VALUE { {{ {{1}}, {{0}}, {{0}}, {{0}}, 0 }} };
};

template <typename TVoidType>
constexpr std::array<OptimalSearchGM<1>, 1> OptimalSearchSchemesGM<0, TVoidType>::VALUE;

template <typename TVoidType>
struct OptimalSearchSchemesGM<1, TVoidType>
{
    static constexpr std::array<OptimalSearchGM<2>, 2> VALUE
    {{
        { {{1, 2}}, {{0, 0}}, {{0, 1}}, {{0, 0}}, 0 },
        { {{2, 1}}, {{0, 1}}, {{0, 1}}, {{0, 0}}, 0 }
    }};
};

template <typename TVoidType>
constexpr std::array<OptimalSearchGM<2>, 2> OptimalSearchSchemesGM<1, TVoidType>::VALUE;

template <typename TVoidType>
struct OptimalSearchSchemesGM<2, TVoidType>
{
    static constexpr std::array<OptimalSearchGM<4>, 3> VALUE
    {{
        { {{1, 2, 3, 4}}, {{0, 0, 1, 1}}, {{0, 0, 2, 2}}, {{0, 0, 0, 0}}, 0 },
        { {{3, 2, 1, 4}}, {{0, 0, 0, 0}}, {{0, 1, 1, 2}}, {{0, 0, 0, 0}}, 0 },
        { {{4, 3, 2, 1}}, {{0, 0, 0, 2}}, {{0, 1, 2, 2}}, {{0, 0, 0, 0}}, 0 }
    }};
};

template <typename TVoidType>
constexpr std::array<OptimalSearchGM<4>, 3> OptimalSearchSchemesGM<2, TVoidType>::VALUE;

template <typename TVoidType>
struct OptimalSearchSchemesGM<3, TVoidType>
{
    static constexpr std::array<OptimalSearchGM<5>, 4> VALUE
    {{
        { {{1, 2, 3, 4, 5}}, {{0, 0, 0, 0, 3}}, {{0, 1, 2, 3, 3}}, {{0, 0, 0, 0, 0}}, 0 },
        { {{2, 3, 4, 5, 1}}, {{0, 0, 0, 2, 2}}, {{0, 1, 2, 2, 3}}, {{0, 0, 0, 0, 0}}, 0 },
        { {{3, 4, 5, 2, 1}}, {{0, 0, 1, 1, 1}}, {{0, 1, 1, 3, 3}}, {{0, 0, 0, 0, 0}}, 0 },
        { {{5, 4, 3, 2, 1}}, {{0, 0, 0, 0, 0}}, {{0, 0, 3, 3, 3}}, {{0, 0, 0, 0, 0}}, 0 }
    }};
};

template <typename TVoidType>
constexpr std::array<OptimalSearchGM<5>, 4> OptimalSearchSchemesGM<3, TVoidType>::VALUE;

template <typename TVoidType>
struct OptimalSearchSchemesGM<4, TVoidType>
{
    static constexpr std::array<OptimalSearchGM<6>, 7> VALUE
    {{
        { {{1, 2, 3, 4, 5, 6}}, {{0, 0, 0, 0, 0, 4}}, {{0, 2, 3, 3, 4, 4}}, {{0, 0, 0, 0, 0, 0}}, 0 },
        { {{3, 4, 5, 6, 2, 1}}, {{0, 0, 0, 1, 4, 4}}, {{0, 0, 1, 1, 4, 4}}, {{0, 0, 0, 0, 0, 0}}, 0 },
        { {{2, 3, 4, 5, 6, 1}}, {{0, 0, 0, 0, 0, 0}}, {{0, 2, 2, 3, 3, 4}}, {{0, 0, 0, 0, 0, 0}}, 0 },
        { {{3, 2, 4, 5, 6, 1}}, {{0, 1, 1, 1, 1, 1}}, {{0, 1, 2, 3, 3, 4}}, {{0, 0, 0, 0, 0, 0}}, 0 },
        { {{4, 3, 2, 5, 6, 1}}, {{0, 0, 2, 2, 2, 2}}, {{0, 0, 2, 3, 3, 4}}, {{0, 0, 0, 0, 0, 0}}, 0 },
        { {{4, 3, 2, 5, 6, 1}}, {{0, 1, 2, 2, 2, 2}}, {{0, 1, 2, 3, 3, 4}}, {{0, 0, 0, 0, 0, 0}}, 0 },
        { {{6, 5, 4, 3, 2, 1}}, {{0, 0, 0, 0, 3, 3}}, {{0, 0, 4, 4, 4, 4}}, {{0, 0, 0, 0, 0, 0}}, 0 }
    }};
};

template <typename TVoidType>
constexpr std::array<OptimalSearchGM<6>, 7> OptimalSearchSchemesGM<4, TVoidType>::VALUE;

// Given the blocklengths (absolute, not cumulative values), assign it to all
// OptimalSearches in a OptimalSearchScheme. The order of blocklength has to be from left to
// right (regarding blocks)
template <size_t nbrBlocks, size_t N>
inline void _optimalSearchSchemeSetBlockLengthGM(std::array<OptimalSearchGM<nbrBlocks>, N> & ss,
                                                 std::vector<uint32_t> const & blocklength)
{
    for (OptimalSearchGM<nbrBlocks> & s : ss)
        for (uint8_t i = 0; i < s.blocklength.size(); ++i)
            s.blocklength[i] = blocklength[s.pi[i]-1] + ((i > 0) ? s.blocklength[i-1] : 0);
}

// requires blocklength to be already set!
template <size_t nbrBlocks, size_t N>
inline void _optimalSearchSchemeInitGM(std::array<OptimalSearchGM<nbrBlocks>, N> & ss)
{
    // check whether 2nd block is on the left or right and choose initialDirection accordingly
    // (more efficient since we do not have to switch directions and thus have better caching performance)
    // for that we need to slightly modify search()
    for (OptimalSearchGM<nbrBlocks> & s : ss)
    {
        s.startPos = 0;
        for (uint8_t i = 0; i < s.pi.size(); ++i)
            if (s.pi[i] < s.pi[0])
                s.startPos += s.blocklength[i] - s.blocklength[i - 1];
    }
}

template <size_t nbrBlocks, size_t N>
inline void _optimalSearchSchemeComputeFixedBlocklengthGM(std::array<OptimalSearchGM<nbrBlocks>, N> & ss, uint32_t const needleLength)
{
    uint8_t blocks = ss[0].pi.size();
    uint32_t blocklength = needleLength / blocks;
    uint8_t rest = needleLength - blocks * blocklength;
    std::vector<uint32_t> blocklengths;
    for (uint8_t i = 0; i < blocks; ++i)
        blocklengths.push_back(blocklength + (i < rest));

    _optimalSearchSchemeSetBlockLengthGM(ss, blocklengths);
    _optimalSearchSchemeInitGM(ss);
}

template <typename TDelegate,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TNeedle,
          size_t nbrBlocks,
          typename TDir>
inline void _optimalSearchSchemeDeletionGM(TDelegate & delegate,
                                           Iter<Index<TText, BidirectionalIndex<TIndex> >,   VSTree<TopDown<TIndexSpec> > > iter,
                                           TNeedle const & needle,
                                           uint32_t const needleLeftPos,
                                           uint32_t const needleRightPos,
                                           uint8_t const errors,
                                           OptimalSearchGM<nbrBlocks> const & s,
                                           uint8_t const blockIndex,
                                           TDir const & /**/)
{
    uint8_t const maxErrorsLeftInBlock = s.u[blockIndex] - errors;
    uint8_t const minErrorsLeftInBlock = (s.l[blockIndex] > errors) ? (s.l[blockIndex] - errors) : 0;

    if (minErrorsLeftInBlock == 0)
    {
        uint8_t const blockIndex2 = std::min(blockIndex + 1, static_cast<uint8_t>(s.u.size()) - 1);
        bool const goToRight2 = s.pi[blockIndex2] > s.pi[blockIndex2 - 1];

        if (goToRight2)
        {
            _optimalSearchSchemeGM(delegate, iter, needle, needleLeftPos, needleRightPos, errors, s, blockIndex2, Rev(),
                                   EditDistance());
        }
        else
        {
            _optimalSearchSchemeGM(delegate, iter, needle, needleLeftPos, needleRightPos, errors, s, blockIndex2, Fwd(),
                                   EditDistance());
        }
    }

    if (maxErrorsLeftInBlock > 0 && goDown(iter, TDir()))
    {
        do
        {
            _optimalSearchSchemeDeletionGM(delegate, iter, needle, needleLeftPos, needleRightPos, errors + 1, s,
                                           blockIndex, TDir());
        } while (goRight(iter, TDir()));
    }
}

template <typename TDelegate,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TNeedle,
          size_t nbrBlocks,
          typename TDir,
          typename TDistanceTag>
inline void _optimalSearchSchemeChildrenGM(TDelegate & delegate,
                                           Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
                                           TNeedle const & needle,
                                           uint32_t const needleLeftPos,
                                           uint32_t const needleRightPos,
                                           uint8_t const errors,
                                           OptimalSearchGM<nbrBlocks> const & s,
                                           uint8_t const blockIndex,
                                           uint8_t const minErrorsLeftInBlock,
                                           TDir const & /**/,
                                           TDistanceTag const & /**/)
{
    constexpr bool goToRight = std::is_same<TDir, Rev>::value;
    if (goDown(iter, TDir()))
    {
        constexpr bool isDna5 = std::is_same<typename Value<TNeedle>::Type, Dna5>::value;
        auto const & needleChar = needle[goToRight ? needleRightPos - 1 : needleLeftPos - 1];

        uint32_t charsLeft = s.blocklength[blockIndex] - (needleRightPos - needleLeftPos - 1);
        do
        {
            bool delta = !ordEqual(parentEdgeLabel(iter, TDir()), needleChar) || (isDna5 && needleChar == Dna5('N'));

            // NOTE (cpockrandt): this might not be optimal yet! we have more edges than in the theoretical model,
            // since we go down an edge before we check whether it can even work out!
            if (!std::is_same<TDistanceTag, EditDistance>::value && minErrorsLeftInBlock > 0 &&
                charsLeft + delta < minErrorsLeftInBlock + 1u) // charsLeft - 1 < minErrorsLeftInBlock - delta
            {
                continue;
            }

            int32_t needleLeftPos2 = needleLeftPos - !goToRight;
            uint32_t needleRightPos2 = needleRightPos + goToRight;

            if (needleRightPos - needleLeftPos == s.blocklength[blockIndex])
            {
                // leave the possibility for one or multiple deletions! therefore, don't change direction, etc!
                SEQAN_IF_CONSTEXPR (std::is_same<TDistanceTag, EditDistance>::value)
                {
                    _optimalSearchSchemeDeletionGM(delegate, iter, needle, needleLeftPos2, needleRightPos2,
                                                   errors + delta, s, blockIndex, TDir());
                }
                else
                {
                    uint8_t blockIndex2 = std::min(blockIndex + 1, static_cast<uint8_t>(s.u.size()) - 1);
                    bool goToRight2 = s.pi[blockIndex2] > s.pi[blockIndex2 - 1];
                    if (goToRight2)
                    {
                        _optimalSearchSchemeGM(delegate, iter, needle, needleLeftPos2, needleRightPos2, errors + delta, s,
                                               blockIndex2, Rev(), TDistanceTag());
                    }
                    else
                    {
                        _optimalSearchSchemeGM(delegate, iter, needle, needleLeftPos2, needleRightPos2, errors + delta, s,
                                               blockIndex2, Fwd(), TDistanceTag());
                    }
                }
            }
            else
            {
                _optimalSearchSchemeGM(delegate, iter, needle, needleLeftPos2, needleRightPos2, errors + delta, s,
                                       blockIndex, TDir(), TDistanceTag());
            }

            // Deletion
            SEQAN_IF_CONSTEXPR (std::is_same<TDistanceTag, EditDistance>::value)
            {
                _optimalSearchSchemeGM(delegate, iter, needle, needleLeftPos, needleRightPos, errors + 1, s, blockIndex,
                                       TDir(), TDistanceTag());
            }
        } while (goRight(iter, TDir()));
    }
}

template <typename TDelegate,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TNeedle,
          size_t nbrBlocks,
          typename TDir,
          typename TDistanceTag>
inline void _optimalSearchSchemeExactGM(TDelegate & delegate,
                                        Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
                                        TNeedle const & needle,
                                        uint32_t const needleLeftPos,
                                        uint32_t const needleRightPos,
                                        uint8_t const errors,
                                        OptimalSearchGM<nbrBlocks> const & s,
                                        uint8_t const blockIndex,
                                        TDir const & /**/,
                                        TDistanceTag const & /**/)
{
    constexpr bool isDna5 = std::is_same<typename Value<TNeedle>::Type, Dna5>::value;
    bool goToRight2 = (blockIndex < s.pi.size() - 1) && s.pi[blockIndex + 1] > s.pi[blockIndex];
    SEQAN_IF_CONSTEXPR (std::is_same<TDir, Rev>::value)
    {
        uint32_t infixPosLeft = needleRightPos - 1;
        uint32_t infixPosRight = needleLeftPos + s.blocklength[blockIndex] - 1;

        // TODO: search first for N and then perform goDown (same for extend and approx in algo.hpp)
        while (infixPosLeft <= infixPosRight)
        {
            if ((isDna5 && needle[infixPosLeft] == Dna5('N')) || !goDown(iter, needle[infixPosLeft], TDir()))
                return;
            ++infixPosLeft;
        }

        if (goToRight2)
        {
            _optimalSearchSchemeGM(delegate, iter, needle, needleLeftPos, infixPosRight + 2, errors, s,
                                   std::min(blockIndex + 1, static_cast<uint8_t>(s.u.size()) - 1), Rev(), TDistanceTag());
        }
        else
        {
            _optimalSearchSchemeGM(delegate, iter, needle, needleLeftPos, infixPosRight + 2, errors, s,
                                   std::min(blockIndex + 1, static_cast<uint8_t>(s.u.size()) - 1), Fwd(), TDistanceTag());
        }
    }
    else
    {
        // has to be signed, otherwise we run into troubles when checking for -1 >= 0u
        int32_t infixPosLeft = needleRightPos - s.blocklength[blockIndex] - 1;
        int32_t infixPosRight = needleLeftPos - 1;

        while (infixPosLeft <= infixPosRight)
        {
            if ((isDna5 && needle[infixPosRight] == Dna5('N')) || !goDown(iter, needle[infixPosRight], TDir()))
                return;
            --infixPosRight;
        }
        if (goToRight2)
        {
            _optimalSearchSchemeGM(delegate, iter, needle, infixPosLeft, needleRightPos, errors, s,
                                   std::min(blockIndex + 1, static_cast<uint8_t>(s.u.size()) - 1), Rev(), TDistanceTag());
        }
        else
        {
            _optimalSearchSchemeGM(delegate, iter, needle, infixPosLeft, needleRightPos, errors, s,
                                   std::min(blockIndex + 1, static_cast<uint8_t>(s.u.size()) - 1), Fwd(), TDistanceTag());
        }
    }
}

template <typename TDelegate,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TNeedle,
          size_t nbrBlocks,
          typename TDir,
          typename TDistanceTag>
inline void _optimalSearchSchemeGM(TDelegate & delegate,
                                   Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > iter,
                                   TNeedle const & needle,
                                   uint32_t const needleLeftPos,
                                   uint32_t const needleRightPos,
                                   uint8_t const errors,
                                   OptimalSearchGM<nbrBlocks> const & s,
                                   uint8_t const blockIndex,
                                   TDir const & /**/,
                                   TDistanceTag const & /**/)
{
    uint8_t const maxErrorsLeftInBlock = s.u[blockIndex] - errors;
    uint8_t const minErrorsLeftInBlock = (s.l[blockIndex] > errors) ? (s.l[blockIndex] - errors) : 0;

    // Done.
    if (minErrorsLeftInBlock == 0 && needleLeftPos == 0 && needleRightPos == length(needle) + 1)
    {
        delegate(iter, needle, errors);
    }
    // Exact search in current block.
    else if (maxErrorsLeftInBlock == 0 && needleRightPos - needleLeftPos - 1 != s.blocklength[blockIndex])
    {
        _optimalSearchSchemeExactGM(delegate, iter, needle, needleLeftPos, needleRightPos, errors, s, blockIndex, TDir(),
                                    TDistanceTag());
    }
    // Approximate search in current block.
    // (s.blocklength[blockIndex]-(needleRightPos-needleLeftPos-(needleLeftPos!=needleRightPos))>=minErrorsLeftInBlock)
    else
    {
        // Insertion
        SEQAN_IF_CONSTEXPR (std::is_same<TDistanceTag, EditDistance>::value)
        {
            constexpr bool goToRight = std::is_same<TDir, Rev>::value;
            int32_t const needleLeftPos2 = needleLeftPos - !goToRight;
            uint32_t const needleRightPos2 = needleRightPos + goToRight;

            if (needleRightPos - needleLeftPos == s.blocklength[blockIndex])
            {
                // leave the possibility for one or multiple deletions! therefore, don't change direction, etc!
                _optimalSearchSchemeDeletionGM(delegate, iter, needle, needleLeftPos2, needleRightPos2, errors + 1, s,
                                               blockIndex, TDir());
            }
            else
            {
                _optimalSearchSchemeGM(delegate, iter, needle, needleLeftPos2, needleRightPos2, errors + 1, s, blockIndex,
                                       TDir(), TDistanceTag());
            }
        }
        _optimalSearchSchemeChildrenGM(delegate, iter, needle, needleLeftPos, needleRightPos, errors, s, blockIndex,
                                       minErrorsLeftInBlock, TDir(), TDistanceTag());
    }
}

template <typename TDelegate,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TNeedle,
          size_t nbrBlocks,
          typename TDistanceTag>
inline void _optimalSearchSchemeGM(TDelegate & delegate,
                                   Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > it,
                                   TNeedle const & needle,
                                   OptimalSearchGM<nbrBlocks> const & s,
                                   TDistanceTag const & /**/)
{
    _optimalSearchSchemeGM(delegate, it, needle, s.startPos, s.startPos + 1, 0, s, 0, Rev(), TDistanceTag());
}

template <typename TDelegate,
          typename TText, typename TIndex, typename TIndexSpec,
          typename TNeedle,
          size_t nbrBlocks, size_t N,
          typename TDistanceTag>
inline void _optimalSearchSchemeGM(TDelegate & delegate,
                                   Iter<Index<TText, BidirectionalIndex<TIndex> >, VSTree<TopDown<TIndexSpec> > > it,
                                   TNeedle const & needle,
                                   std::array<OptimalSearchGM<nbrBlocks>, N> const & ss,
                                   TDistanceTag const & /**/)
{
    for (auto & s : ss)
        _optimalSearchSchemeGM(delegate, it, needle, s, TDistanceTag());
}

}
