#pragma once

#include <time.h>
#include <sys/time.h>

#include <seqan/index.h>

using namespace seqan;

inline auto retrieveDirectoryInformationLine(CharString const & info)
{
    std::string const row = toCString(info);
    auto const firstSeparator = row.find(';', 0);
    auto const secondSeparator = row.find(';', firstSeparator + 1);
    std::string const fastaFile = row.substr(0, firstSeparator);
    uint64_t const length = std::stoi(row.substr(firstSeparator + 1, secondSeparator - firstSeparator - 1));
    std::string const chromName = row.substr(secondSeparator + 1);
    return std::make_tuple(fastaFile, length, chromName);
}

template <typename TResult, typename TPosition, typename TLimits>
inline void myPosLocalize(TResult & result, TPosition const & pos, TLimits const & limits) {
    typedef typename Iterator<TLimits const, Standard>::Type TIter;
    TIter _begin = begin(limits, Standard());
    TIter _upper = std::upper_bound(_begin, end(limits, Standard()), pos) - 1;
    result.i1 = difference(_begin, _upper);
    result.i2 = pos - *_upper;
}

inline double get_wall_time()
{
    struct timeval time;
    if (gettimeofday(&time, NULL))
        return 0;
    return static_cast<double>(time.tv_sec) + static_cast<double>(time.tv_usec) * .000001;
}

template <typename TSpec = void, typename TLengthSum = size_t, unsigned LEVELS = 2, unsigned WORDS_PER_BLOCK = 1>
struct GemMapFastFMIndexConfig
{
    typedef TLengthSum                                                                          LengthSum;
    typedef Levels<TSpec, LevelsPrefixRDConfig<LengthSum, Alloc<>, LEVELS, WORDS_PER_BLOCK> >   Bwt;
    typedef Levels<TSpec, LevelsRDConfig<LengthSum, Alloc<>, LEVELS, WORDS_PER_BLOCK> >         Sentinels;

    static unsigned SAMPLING;
};

template <typename TLengthSum>
using TGemMapFastFMIndexConfig = GemMapFastFMIndexConfig<void, TLengthSum, 2, 1>;

template <typename TFMIndexConfig>
using TBiIndexConfig = BidirectionalIndex<FMIndex<void, TFMIndexConfig> >;

namespace seqan {

template <typename TSeqNo, typename TSeqPos>
struct SizeSpec_ {};

template <typename TString, typename TSeqNo, typename TSeqPos>
struct SAValue<StringSet<TString, Owner<ConcatDirect<SizeSpec_<TSeqNo, TSeqPos> > > > >
{
    typedef Pair<TSeqNo, TSeqPos, Pack> Type;
};

} // namespace seqan

struct SearchParams
{
    unsigned length;
    unsigned overlap;
    unsigned threads;
    // bool indels;
    bool revCompl;
    bool excludePseudo;
};

std::string mytime()
{
    auto r = time(nullptr);
    auto c = ctime(&r);
    std::string buf(c);
    buf.insert(0, "[");
    buf.append("] ");
    buf.erase(remove(buf.begin(), buf.end(), '\n'), buf.end());
    return buf;
}

template <bool outputProgress>
inline void printProgress(uint64_t &, uint64_t const, uint64_t const);

template <>
inline void printProgress<false>(uint64_t &, uint64_t const, uint64_t const)
{ }

template <>
inline void printProgress<true>(uint64_t & progress_count, uint64_t const progress_step, uint64_t const progress_max)
{
    #pragma omp atomic
    ++progress_count;
    if (omp_get_thread_num() == 0 && (progress_count & progress_step) == 0)
    {
        float progress = static_cast<float>(progress_count)/progress_max;
        std::cout << "\rProgress: " << (truncf(progress*10000)/100) << "%   " << std::flush;
    }
}

template <bool outputProgress>
inline void initProgress(uint64_t &, uint64_t &, uint64_t &, uint64_t const, uint64_t const);

template <>
inline void initProgress<false>(uint64_t &, uint64_t &, uint64_t &, uint64_t const, uint64_t const)
{ }

template <>
inline void initProgress<true>(uint64_t & progressCount, uint64_t & progressStep, uint64_t & progressMax,
                               uint64_t const stepSize, uint64_t const numberOfKmers)
{
    progressCount = 0;
    progressMax = (numberOfKmers + stepSize - 1) / stepSize; // = ceil(numberOfKmers / stepSize), i.e. loop iterations
    progressStep = 511; // Print every 512 loop iterators.
}
