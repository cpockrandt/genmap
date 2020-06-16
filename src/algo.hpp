#include "find2_index_approx.hpp"

using namespace seqan;

template <typename TChar>
using ModComplement = ModView<FunctorComplement<TChar> >;
template <typename TText>
using ModRevCompl = ModifiedString<ModifiedString<TText, ModComplement<typename Value<TText>::Type>>, ModReverse>;

template <typename TMappVector, typename TCumChromLengths>
void resetLimits(TMappVector & c, unsigned const kmerLength, TCumChromLengths const & cumChromLengths)
{
    // skip first, since the first cumulative length is 0
    for (uint64_t i = 1; i < length(cumChromLengths); ++i)
    {
        // std::min(): make sure that we don't mess up if a sequence is shorter than K
        for (uint64_t j = 1; j < std::min<uint64_t>(kmerLength, cumChromLengths[i] - cumChromLengths[i - 1] + 1); ++j)
        {
            c[cumChromLengths[i] - j] = 0;
        }
    }
}

// TODO: avoid signed integers

template <bool reportExactMatch, unsigned maxErrors, typename TBiIter, typename TValue, typename TText>
inline void extendExact(TBiIter it, std::vector<TValue> & hits, std::vector<typename TBiIter::TFwdIndexIter> & itExact,
                        std::vector<std::vector<typename TBiIter::TFwdIndexIter> > & itAll,
                        TText const & text, unsigned const length,
                        uint64_t a, uint64_t b, // searched interval
                        uint64_t ab, uint64_t bb, // entire interval
                        bool const csvComputation)
{
    constexpr bool isDna5 = std::is_same<typename Value<TText>::Type, Dna5>::value;

    constexpr uint64_t max_val = std::numeric_limits<TValue>::max();

    if (b - a + 1 == length)
    {
        SEQAN_IF_CONSTEXPR (reportExactMatch && maxErrors == 0)
        {
            itExact[a-ab] = it.fwdIter;
        }
        if (csvComputation)
        {
            itAll[a-ab].push_back(it.fwdIter);
        }
        hits[a-ab] = std::min((uint64_t) countOccurrences(it) + hits[a-ab], max_val);
        return;
    }
    //if (b + 1 <= bb)
    //{
        TBiIter it2 = it;
        uint64_t brm = a + length - 1;
        uint64_t b_new = b + (((brm - b) + 2 - 1) >> 1); // ceil((bb - b)/2)
        if (b_new <= bb)
        {
            bool success = true;
            for (uint64_t i = b + 1; i <= b_new && success; ++i)
            {
                success = (!isDna5 || text[i] != Dna5('N')) && goDown(it2, text[i], Rev());
            }
            if (success)
                extendExact<reportExactMatch, maxErrors>(it2, hits, itExact, itAll, text, length, a, b_new, ab, bb, csvComputation);
        }
    //}

    if (a - 1 >= ab)
    {
        int64_t alm = b + 1 - length;
        uint64_t a_new = alm + std::max<int64_t>(((a - alm) - 1) >> 1, 0);
        for (int64_t i = a - 1; i >= static_cast<int64_t>(a_new); --i)
        {
            if((isDna5 && text[i] == Dna5('N')) || !goDown(it, text[i], Fwd()))
                return;
        }
        extendExact<reportExactMatch, maxErrors>(it, hits, itExact, itAll, text, length, a_new, b, ab, bb, csvComputation);
    }
}

// forward
template <bool reportExactMatch, unsigned maxErrors, typename TBiIter, typename TValue, typename TText>
inline void extend(TBiIter it, std::vector<TValue> & hits, std::vector<typename TBiIter::TFwdIndexIter> & itExact,
                   std::vector<std::vector<typename TBiIter::TFwdIndexIter> > & itAll,
                   unsigned errorsLeft, TText const & text, unsigned const length,
                   uint64_t a, uint64_t b, // searched interval
                   uint64_t ab, uint64_t bb, // entire interval
                   bool const csvComputation);

template <bool reportExactMatch, unsigned maxErrors, typename TBiIter, typename TValue, typename TText>
inline void approxSearch(TBiIter it, std::vector<TValue> & hits, std::vector<typename TBiIter::TFwdIndexIter> & itExact,
                         std::vector<std::vector<typename TBiIter::TFwdIndexIter> > & itAll,
                         unsigned errorsLeft, TText const & text, unsigned const length,
                         uint64_t a, uint64_t b, // searched interval
                         uint64_t ab, uint64_t bb, // entire interval
                         uint64_t b_new,
                         Rev const &, bool const csvComputation)
{
    constexpr bool isDna5 = std::is_same<typename Value<TText>::Type, Dna5>::value;

    if (b == b_new)
    {
        extend<reportExactMatch, maxErrors>(it, hits, itExact, itAll, errorsLeft, text, length, a, b, ab, bb, csvComputation);
        return;
    }
    if (errorsLeft > 0)
    {
        if (goDown(it, Rev()))
        {
            do {
                bool delta = !ordEqual(parentEdgeLabel(it, Rev()), text[b + 1])
                             || (isDna5 && text[b + 1] == Dna5('N'));
                approxSearch<reportExactMatch, maxErrors>(it, hits, itExact, itAll, errorsLeft - delta, text, length, a, b + 1, ab, bb, b_new, Rev(), csvComputation);
            } while (goRight(it, Rev()));
        }
    }
    else
    {
        for (uint64_t i = b + 1; i <= b_new; ++i)
        {
            if ((isDna5 && text[i] == Dna5('N')) || !goDown(it, text[i], Rev()))
                return;
        }
        extendExact<reportExactMatch, maxErrors>(it, hits, itExact, itAll, text, length, a, b_new, ab, bb, csvComputation);
    }
}
template <bool reportExactMatch, unsigned maxErrors, typename TBiIter, typename TValue, typename TText>
inline void approxSearch(TBiIter it, std::vector<TValue> & hits, std::vector<typename TBiIter::TFwdIndexIter> & itExact,
                         std::vector<std::vector<typename TBiIter::TFwdIndexIter> > & itAll,
                         unsigned errorsLeft, TText const & text, unsigned const length,
                         uint64_t a, uint64_t b, // searched interval
                         uint64_t ab, uint64_t bb, // entire interval
                         uint64_t a_new,
                         Fwd const & /*tag*/, bool const csvComputation)
{
    constexpr bool isDna5 = std::is_same<typename Value<TText>::Type, Dna5>::value;

    if (a == a_new)
    {
        extend<reportExactMatch, maxErrors>(it, hits, itExact, itAll, errorsLeft, text, length, a, b, ab, bb, csvComputation);
        return;
    }
    if (errorsLeft > 0)
    {
        if (goDown(it, Fwd()))
        {
            do {
                bool delta = !ordEqual(parentEdgeLabel(it, Fwd()), text[a - 1])
                             || (isDna5 && text[a - 1] == Dna5('N'));
                approxSearch<reportExactMatch, maxErrors>(it, hits, itExact, itAll, errorsLeft - delta, text, length, a - 1, b, ab, bb, a_new, Fwd(), csvComputation);
            } while (goRight(it, Fwd()));
        }
    }
    else
    {
        for (int64_t i = a - 1; i >= static_cast<int64_t>(a_new); --i)
        {
            if ((isDna5 && text[i] == Dna5('N')) || !goDown(it, text[i], Fwd()))
                return;
        }
        extendExact<reportExactMatch, maxErrors>(it, hits, itExact, itAll, text, length, a_new, b, ab, bb, csvComputation);
    }
}

template <bool reportExactMatch, unsigned maxErrors, typename TBiIter, typename TValue, typename TText>
inline void extend(TBiIter it, std::vector<TValue> & hits, std::vector<typename TBiIter::TFwdIndexIter> & itExact,
                   std::vector<std::vector<typename TBiIter::TFwdIndexIter> > & itAll,
                   unsigned errorsLeft, TText const & text, unsigned const length,
                   uint64_t a, uint64_t b, // searched interval
                   uint64_t ab, uint64_t bb, // entire interval
                   bool const csvComputation)
{
    constexpr uint64_t max_val = std::numeric_limits<TValue>::max();

    if (errorsLeft == 0)
    {
        extendExact<reportExactMatch, maxErrors>(it, hits, itExact, itAll, text, length, a, b, ab, bb, csvComputation);
        return;
    }
    if (b - a + 1 == length)
    {
        SEQAN_IF_CONSTEXPR (reportExactMatch)
        {
            if (maxErrors == errorsLeft)
                itExact[a-ab] = it.fwdIter;
        }
        if (csvComputation)
        {
            itAll[a-ab].push_back(it.fwdIter);
        }
        hits[a-ab] = std::min((uint64_t) countOccurrences(it) + hits[a-ab], max_val);
        return;
    }
    //if (b + 1 <= bb)
    //{
        uint64_t brm = a + length - 1;
        uint64_t b_new = b + (((brm - b) + 2 - 1) >> 1); // ceil((bb - b)/2)
        if (b_new <= bb)
        {
            approxSearch<reportExactMatch, maxErrors>(it, hits, itExact, itAll, errorsLeft, text, length,
                         a, b, // searched interval
                         ab, bb, // entire interval
                         b_new,
                         Rev(), csvComputation);
        }
    //}

    if (a - 1 >= ab)
    {
        int64_t alm = b + 1 - length;
        uint64_t a_new = alm + std::max<int64_t>(((a - alm) - 1) >> 1, 0);
        approxSearch<reportExactMatch, maxErrors>(it, hits, itExact, itAll, errorsLeft, text, length,
                     a, b, // searched interval
                     ab, bb, // entire interval
                     a_new,
                     Fwd(), csvComputation);
    }
}

// computes a block of adjacent k-mers at once
template <unsigned errors, typename TIndex, typename TText, typename TContainer, typename TChromosomeLengths, typename TLocations, typename TMapping, typename TLimits>
inline void computeMappabilitySingleBlock(TIndex & index, TText const & text, TContainer & c, SearchParams const & params,
                                          bool const directory, TChromosomeLengths const & chromLengths, TLocations & locations, TMapping const & mappingSeqIdFile,
                                          uint64_t const i, uint64_t const j, uint64_t const textLength, TChromosomeLengths const & chromCumLengths, TLimits const & limits,
                                          std::vector<std::pair<uint64_t, uint64_t>> const & intervals, unsigned const overlap, bool const completeSameKmers, bool const csvComputation)
{
    typedef typename TContainer::value_type TValue;
    typedef Iter<TIndex, VSTree<TopDown<> > > TBiIter;

    // overlap is the length of the infix!
    uint64_t maxPos = std::min(i + params.length - overlap, textLength - params.length) + 1;
    if (maxPos > j)
        maxPos = j;

    // Skip leading and trailing precomputed k-mer frequencies
    uint64_t beginPos = i;
    while (beginPos < maxPos && c[beginPos] != 0)
        ++beginPos;

    uint64_t endPos = maxPos; // endPos is excluding, i.e. [beginPos, endPos)
    while (i > 0 && endPos - 1 >= i && c[endPos - 1] != 0) // we do not check for i == 0 to avoid an underflow.
        --endPos;

    if (beginPos < endPos)
    {
        uint64_t overlap = params.length - (endPos - beginPos) + 1;

        auto scheme = OptimalSearchSchemesGM<errors>::VALUE;
        _optimalSearchSchemeComputeFixedBlocklengthGM(scheme, overlap);

        std::vector<typename TBiIter::TFwdIndexIter> itExact(endPos - beginPos);
        std::vector<TValue> hits(endPos - beginPos, 0);
        std::vector<std::vector<typename TBiIter::TFwdIndexIter> > itAll(endPos - beginPos);
        std::vector<std::vector<typename TBiIter::TFwdIndexIter> > itAllrevCompl(endPos - beginPos);

        auto const & needles = infix(text, beginPos, beginPos + params.length + (endPos - beginPos) - 1);
        auto const & needlesOverlap = infix(text, beginPos + params.length - overlap, beginPos + params.length);
        using TNeedlesOverlap = decltype(needlesOverlap);

        uint64_t const bb = std::min(textLength - 1, params.length - 1 + params.length - overlap);

        auto delegate = [&hits, &itExact, &itAll, bb, overlap, &params, &needles, csvComputation](
            TBiIter it, TNeedlesOverlap const & /*read*/, unsigned const errors_spent)
        {
            // TODO: we could turn reporting of exact iterators off at compile time by setting reportExactMatch = false if opt.directory is true. Evaluate binary size vs. performance.
            // WARNING: if it is computed on the directory, csvComputation currently still needs the exact matches (can be updated down below)

            if (errors_spent == 0)
            {
                extend<true, errors>(it, hits, itExact, itAll, errors - errors_spent, needles, params.length,
                    params.length - overlap, params.length - 1, // searched interval
                    0, bb, // entire interval
                    csvComputation);
            }
            else
            {
                extend<false, errors>(it, hits, itExact, itAll, errors - errors_spent, needles, params.length,
                    params.length - overlap, params.length - 1, // searched interval
                    0, bb, // entire interval
                    csvComputation);
            }
        };

        if (params.revCompl)
        {
            ModRevCompl<typename std::remove_reference<decltype(needles)>::type> needlesRevCompl(needles);
            ModRevCompl<typename std::remove_reference<decltype(needlesOverlap)>::type> needlesRevComplOverlap(needlesOverlap);
            using TNeedlesRevComplOverlap = decltype(needlesRevComplOverlap);

            // TODO: could store the exact hits as well and use these values!
            auto delegateRevCompl = [&hits, &itExact, &itAllrevCompl, bb, overlap, &params, &needlesRevCompl, csvComputation](
                TBiIter it, TNeedlesRevComplOverlap const & /*read*/, unsigned const errors_spent)
            {
                extend<false, errors>(it, hits, itExact, itAllrevCompl, errors - errors_spent, needlesRevCompl, params.length,
                    params.length - overlap, params.length - 1, // searched interval
                    0, bb, // entire interval
                    csvComputation);
            };

            TBiIter it(index);
            _optimalSearchSchemeGM(delegateRevCompl, it, needlesRevComplOverlap, scheme, HammingDistance());

            // hits of the reverse-complement are stored in reversed order.
            std::reverse(hits.begin(), hits.end());
        }

        TBiIter it(index);
        _optimalSearchSchemeGM(delegate, it, needlesOverlap, scheme, HammingDistance());
        for (uint64_t j = beginPos; j < endPos; ++j)
        {
            if (csvComputation)
            {
                using TLocation = typename TLocations::key_type;
                using TEntry = std::pair<TLocation, std::pair<std::vector<TLocation>, std::vector<TLocation> > >;

                TEntry entry;

                uint64_t size = 0;
                for (auto const & iterator : itAll[j - beginPos])
                    size += countOccurrences(iterator);
                entry.second.first.reserve(size);

                size = 0;
                for (auto const & iterator : itAllrevCompl[j - beginPos])
                    size += countOccurrences(iterator);
                entry.second.second.reserve(size);

                for (auto const & iterator : itAll[j - beginPos])
                {
                    for (auto const & occ : getOccurrences(iterator))
                    {
                        entry.second.first.push_back(occ);
                    }
                }
                // sorting is needed for output when multiple fasta files are indexed and the locations need to be separated by filename.
                std::sort(entry.second.first.begin(), entry.second.first.end());

                // NOTE: vector has to be iterated over in reverse order (compared to itAll)
                // for (auto const & iterator : itAllrevCompl[j - beginPos])
                for (auto const & iterator : itAllrevCompl[endPos - 1 - j])
                {
                    for (auto const & occ : getOccurrences(iterator))
                    {
                        entry.second.second.push_back(occ);
                    }
                }
                // sorting is needed for output when multiple fasta files are indexed and the locations need to be separated by filename.
                std::sort(entry.second.second.begin(), entry.second.second.end());

                // overwrite frequency vector
                if (params.excludePseudo)
                {
                    std::set<typename Value<TLocation, 1>::Type> distinct_sequences;
                    for (auto const & location : entry.second.first) // forward strand
                        distinct_sequences.emplace(mappingSeqIdFile[location.i1]);
                    assert(entry.second.second.size() == 0 || params.revCompl);
                    for (auto const & location : entry.second.second) // reverse strand
                        distinct_sequences.emplace(mappingSeqIdFile[location.i1]);

                    hits[j - beginPos] = distinct_sequences.size();

                    // NOTE: If you want to filter certain k-mers in the csv file based on the mappability value
                    // (with respect to --exclude-pseudo) you can unset 'entry' here.
                }

                if (!directory && countOccurrences(itExact[j - beginPos]) > 1)
                {
                    for (auto const & exact_occ : getOccurrences(itExact[j - beginPos]))
                    {
                        if (static_cast<int64_t>(exact_occ.i2) <= static_cast<int64_t>(chromLengths[exact_occ.i1]) - params.length)
                        {
                            #pragma omp critical
                            locations.emplace(exact_occ, entry.second);
                        }
                    }
                }
                // is there at least a hit on the forward or the reverse strand? This is needed for Dna5
                else if (entry.second.first.size() + entry.second.second.size() > 0)
                {
                    myPosLocalize(entry.first, j, chromCumLengths); // TODO: inefficient for read data sets
                    if (static_cast<int64_t>(entry.first.i2) <= static_cast<int64_t>(chromLengths[entry.first.i1]) - params.length)
                    {
                        #pragma omp critical
                        locations.emplace(entry);
                    }
                }
            }

            if (!directory && (intervals.empty() || completeSameKmers) && countOccurrences(itExact[j - beginPos]) > 1) // guaranteed to exist, since there has to be at least one match!
            {
                for (auto const & occ : getOccurrences(itExact[j-beginPos]))
                {
                    auto const occ_pos = posGlobalize(occ, limits);
                    c[occ_pos] = hits[j - beginPos];
                }
            }
            else
            {
                c[j] = hits[j - beginPos];
            }
        }
    }
}

template <unsigned errors, typename TIndex, typename TText, typename TContainer, typename TChromosomeLengths, typename TLocations, typename TMapping>
inline void computeMappability(TIndex & index, TText const & text, TContainer & c, SearchParams const & params,
                               bool const directory, TChromosomeLengths const & chromLengths, TChromosomeLengths const & chromCumLengths, TLocations & locations,
                               TMapping const & mappingSeqIdFile, std::vector<std::pair<uint64_t, uint64_t>> const & intervals,
                               bool & completeSameKmers,
                               uint64_t const currentFileNo, uint64_t const totalFileNo, bool const csvComputation)
{
    auto const & limits = stringSetLimits(indexText(index));
    uint64_t const textLength = length(text);
    uint64_t const numberOfKmers = textLength - params.length + 1;
    uint64_t const overlap = params.overlap;
    uint64_t const stepSize = params.length - overlap + 1; // Number of overlapping k-mers searched at once

    completeSameKmers = false;

    if (intervals.empty())
    {
        // Number of loop iterations assigned to a thread at once
        // It should be significantly smaller than the number of loop iterations (numberOfKmers/stepSize), since
        // the running time of different loop iterations can vary vastly (e.g., repeats are slower than unique regions).
        // This leads to an unused variable warning in Clang
        #pragma clang diagnostic push
        #pragma clang diagnostic ignored "-Wunused"
        uint64_t const chunkSize = std::max<uint64_t>(1, numberOfKmers / (stepSize * params.threads * 50));
        #pragma clang diagnostic pop

        uint64_t progressCount, progressMax, progressStep;
        initProgress<outputProgress>(progressCount, progressStep, progressMax, stepSize, numberOfKmers);

        #pragma omp parallel for schedule(dynamic, chunkSize) num_threads(params.threads)
        for (uint64_t i = 0; i < numberOfKmers; i += stepSize)
        {
            computeMappabilitySingleBlock<errors>(index, text, c, params, directory, chromLengths, locations, mappingSeqIdFile, i, i + stepSize, textLength, chromCumLengths, limits, intervals, overlap, true, csvComputation);
            printProgress<outputProgress>(progressCount, progressStep, progressMax, currentFileNo, totalFileNo);
        }
    }
    else
    {
        uint64_t interval_sum = 0;
        std::vector<std::pair<uint64_t, uint64_t>> intervals_details;
        for (auto interval = intervals.begin(); interval < intervals.end(); ++interval)
        {
            interval_sum += (*interval).second - (*interval).first;
            for (uint64_t i = (*interval).first; i < (*interval).second; i += stepSize)
            {
                intervals_details.emplace_back(std::make_pair(i, std::min(i + stepSize, (*interval).second)));
            }
        }

        // when bed with subset for mappability is provided with at least 50% of the genome selected
        // use optimization in algorithm (copy mappability of same k-mers)
        float const fraction = static_cast<float>(interval_sum) / textLength;
        completeSameKmers = fraction > 0.5f;

        // This leads to an unused variable warning in Clang
        #pragma clang diagnostic push
        #pragma clang diagnostic ignored "-Wunused"
        uint64_t const chunkSize = std::max<uint64_t>(1, intervals_details.size() / (params.threads * 50));
        #pragma clang diagnostic pop

        uint64_t progressCount, progressMax, progressStep;
        initProgress<outputProgress>(progressCount, progressStep, progressMax, 1, intervals_details.size());

        // NOTE: chunksize for scheduling would depend on number of intervals, size of intervals, deviation of interval sizes, etc.
        // Hence, for simplicity we do not suggest a chunk size
        #pragma omp parallel for schedule(dynamic, chunkSize) num_threads(params.threads)
        for (auto interval = intervals_details.begin(); interval < intervals_details.end(); ++interval)
        {
            computeMappabilitySingleBlock<errors>(index, text, c, params, directory, chromLengths, locations, mappingSeqIdFile, (*interval).first, (*interval).second, textLength, chromCumLengths, limits, intervals, overlap, completeSameKmers, csvComputation);
            printProgress<outputProgress>(progressCount, progressStep, progressMax, currentFileNo, totalFileNo);
        }
    }

    // The algorithm searches k-mers in the concatenation of all strings in the fasta file (e.g. chromosomes).
    // Hence, it also searches k-mers that overlap two strings that actually do not exist.
    // At the end we overwrite the frequency of those k-mers with 0.
    // TODO: k-mers spanning two strings should not be searched if there are many short strings (i.e., fasta of reads).
    resetLimits(c, params.length, chromCumLengths);
}
