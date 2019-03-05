#include "find2_index_approx.hpp"

using namespace seqan;

// template <typename TMappVector, typename TChromosomeLength>
// void resetLimits(TMappVector const &, unsigned const, TChromosomeLength const)
// { }

template <bool csvComputation, typename TMappVector, typename TChromLengths, typename TCumChromLengths, typename TLocations>
void resetLimits(TMappVector & c, unsigned const kmerLength, TChromLengths const & chromLengths, TCumChromLengths const & cumChromLengths, TLocations & locations)
{
    using TLocation = typename TLocations::key_type;
    using TEntry = std::pair<TLocation, std::pair<std::vector<TLocation>, std::vector<TLocation> > >;

    // skip first, since the first cumulative length is 0
    for (uint64_t i = 1; i < length(cumChromLengths); ++i)
    {
        for (uint64_t j = 1; j < kmerLength; ++j)
        {
            c[cumChromLengths[i] - j] = 0;
        }

        // TODO: if sequences are shorter that the kmer (and possibly span more than 2 sequences), this will lead to errors!
        // Remove csv entries for kmers overlapping multiple sequences accidentally.
        SEQAN_IF_CONSTEXPR (csvComputation)
        {
            // for single fasta files it is guaranteed that every entry exists and it only has to be resetted if the kmer is overlapping
            // for directories some kmers might be missing

            // add empty entries in csv for sequences that are shorter that K and reset the last k-1 entries of each sequence
            TEntry entry;
            entry.first.i1 = i - 1;
            entry.first.i2 = (chromLengths[i - 1] >= kmerLength) ? (chromLengths[i - 1] - kmerLength + 1) : 0;
            while (entry.first.i2 < chromLengths[i - 1])
            {
                auto const insertPos = locations.lower_bound(entry.first);
                if (insertPos != locations.end() && !(locations.key_comp()(entry.first, insertPos->first)))
                {
                    insertPos->second.first.clear();
                    insertPos->second.second.clear();
                }
                else
                {
                    locations.insert(insertPos, entry);
                }
                ++entry.first.i2;
            }
        }
    }
}

// TODO: avoid signed integers

template <bool reportExactMatch, bool csvComputation, unsigned maxErrors, typename TBiIter, typename TValue, typename TText>
inline void extendExact(TBiIter it, std::vector<TValue> & hits, std::vector<typename TBiIter::TFwdIndexIter> & itExact,
                        std::vector<std::vector<typename TBiIter::TFwdIndexIter> > & itAll,
                        TText const & text, unsigned const length,
                        uint64_t a, uint64_t b, // searched interval
                        uint64_t ab, uint64_t bb) // entire interval
{
    constexpr bool isDna5 = std::is_same<typename Value<TText>::Type, Dna5>::value;

    constexpr uint64_t max_val = std::numeric_limits<TValue>::max();

    if (b - a + 1 == length)
    {
        SEQAN_IF_CONSTEXPR (reportExactMatch && maxErrors == 0)
        {
            itExact[a-ab] = it.fwdIter;
        }
        SEQAN_IF_CONSTEXPR (csvComputation)
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
                extendExact<reportExactMatch, csvComputation, maxErrors>(it2, hits, itExact, itAll, text, length, a, b_new, ab, bb);
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
        extendExact<reportExactMatch, csvComputation, maxErrors>(it, hits, itExact, itAll, text, length, a_new, b, ab, bb);
    }
}

// forward
template <bool reportExactMatch, bool csvComputation, unsigned maxErrors, typename TBiIter, typename TValue, typename TText>
inline void extend(TBiIter it, std::vector<TValue> & hits, std::vector<typename TBiIter::TFwdIndexIter> & itExact,
                   std::vector<std::vector<typename TBiIter::TFwdIndexIter> > & itAll,
                   unsigned errorsLeft, TText const & text, unsigned const length,
                   uint64_t a, uint64_t b, // searched interval
                   uint64_t ab, uint64_t bb); // entire interval

template <bool reportExactMatch, bool csvComputation, unsigned maxErrors, typename TBiIter, typename TValue, typename TText>
inline void approxSearch(TBiIter it, std::vector<TValue> & hits, std::vector<typename TBiIter::TFwdIndexIter> & itExact,
                         std::vector<std::vector<typename TBiIter::TFwdIndexIter> > & itAll,
                         unsigned errorsLeft, TText const & text, unsigned const length,
                         uint64_t a, uint64_t b, // searched interval
                         uint64_t ab, uint64_t bb, // entire interval
                         uint64_t b_new,
                         Rev const &)
{
    constexpr bool isDna5 = std::is_same<typename Value<TText>::Type, Dna5>::value;

    if (b == b_new)
    {
        extend<reportExactMatch, csvComputation, maxErrors>(it, hits, itExact, itAll, errorsLeft, text, length, a, b, ab, bb);
        return;
    }
    if (errorsLeft > 0)
    {
        if (goDown(it, Rev()))
        {
            do {
                bool delta = !ordEqual(parentEdgeLabel(it, Rev()), text[b + 1])
                             || (isDna5 && text[b + 1] == Dna5('N'));
                approxSearch<reportExactMatch, csvComputation, maxErrors>(it, hits, itExact, itAll, errorsLeft - delta, text, length, a, b + 1, ab, bb, b_new, Rev());
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
        extendExact<reportExactMatch, csvComputation, maxErrors>(it, hits, itExact, itAll, text, length, a, b_new, ab, bb);
    }
}
template <bool reportExactMatch, bool csvComputation, unsigned maxErrors, typename TBiIter, typename TValue, typename TText>
inline void approxSearch(TBiIter it, std::vector<TValue> & hits, std::vector<typename TBiIter::TFwdIndexIter> & itExact,
                         std::vector<std::vector<typename TBiIter::TFwdIndexIter> > & itAll,
                         unsigned errorsLeft, TText const & text, unsigned const length,
                         uint64_t a, uint64_t b, // searched interval
                         uint64_t ab, uint64_t bb, // entire interval
                         uint64_t a_new,
                         Fwd const & /*tag*/)
{
    constexpr bool isDna5 = std::is_same<typename Value<TText>::Type, Dna5>::value;

    if (a == a_new)
    {
        extend<reportExactMatch, csvComputation, maxErrors>(it, hits, itExact, itAll, errorsLeft, text, length, a, b, ab, bb);
        return;
    }
    if (errorsLeft > 0)
    {
        if (goDown(it, Fwd()))
        {
            do {
                bool delta = !ordEqual(parentEdgeLabel(it, Fwd()), text[a - 1])
                             || (isDna5 && text[a - 1] == Dna5('N'));
                approxSearch<reportExactMatch, csvComputation, maxErrors>(it, hits, itExact, itAll, errorsLeft - delta, text, length, a - 1, b, ab, bb, a_new, Fwd());
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
        extendExact<reportExactMatch, csvComputation, maxErrors>(it, hits, itExact, itAll, text, length, a_new, b, ab, bb);
    }
}

template <bool reportExactMatch, bool csvComputation, unsigned maxErrors, typename TBiIter, typename TValue, typename TText>
inline void extend(TBiIter it, std::vector<TValue> & hits, std::vector<typename TBiIter::TFwdIndexIter> & itExact,
                   std::vector<std::vector<typename TBiIter::TFwdIndexIter> > & itAll,
                   unsigned errorsLeft, TText const & text, unsigned const length,
                   uint64_t a, uint64_t b, // searched interval
                   uint64_t ab, uint64_t bb) // entire interval
{
    constexpr uint64_t max_val = std::numeric_limits<TValue>::max();

    if (errorsLeft == 0)
    {
        extendExact<reportExactMatch, csvComputation, maxErrors>(it, hits, itExact, itAll, text, length, a, b, ab, bb);
        return;
    }
    if (b - a + 1 == length)
    {
        SEQAN_IF_CONSTEXPR (reportExactMatch)
        {
            if (maxErrors == errorsLeft)
                itExact[a-ab] = it.fwdIter;
        }
        SEQAN_IF_CONSTEXPR (csvComputation)
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
            approxSearch<reportExactMatch, csvComputation, maxErrors>(it, hits, itExact, itAll, errorsLeft, text, length,
                         a, b, // searched interval
                         ab, bb, // entire interval
                         b_new,
                         Rev()
            );
        }
    //}

    if (a - 1 >= ab)
    {
        int64_t alm = b + 1 - length;
        uint64_t a_new = alm + std::max<int64_t>(((a - alm) - 1) >> 1, 0);
        approxSearch<reportExactMatch, csvComputation, maxErrors>(it, hits, itExact, itAll, errorsLeft, text, length,
                     a, b, // searched interval
                     ab, bb, // entire interval
                     a_new,
                     Fwd()
        );
    }
}

template <unsigned errors, bool csvComputation, typename TIndex, typename TText, typename TContainer, typename TChromosomeLengths, typename TLocations, typename TMapping>
inline void computeMappability(TIndex & index, TText const & text, TContainer & c, SearchParams const & params,
                               bool const directory, TChromosomeLengths const & chromLengths, TLocations & locations, TMapping const & mappingSeqIdFile)
{
    typedef typename TContainer::value_type TValue;
    typedef Iter<TIndex, VSTree<TopDown<> > > TBiIter;

    TChromosomeLengths chromCumLengths;
    {
        uint64_t _cumLength = 0;
        appendValue(chromCumLengths, 0);
        for (uint64_t i = 0; i < length(chromLengths); ++i)
        {
            _cumLength += chromLengths[i];
            appendValue(chromCumLengths, _cumLength);
        }
    }

    auto const & limits = stringSetLimits(indexText(index));
    uint64_t const textLength = length(text);
    uint64_t const numberOfKmers = textLength - params.length + 1;
    uint64_t const stepSize = params.length - params.overlap + 1; // Number of overlapping k-mers searched at once

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
        // overlap is the length of the infix!
        uint64_t maxPos = std::min(i + params.length - params.overlap, textLength - params.length) + 1;

        // Skip leading and trailing precomputed k-mer frequencies
        uint64_t beginPos = i;
        while (beginPos < maxPos && c[beginPos] != 0)
            ++beginPos;

        uint64_t endPos = maxPos; // endPos is excluding, i.e. [beginPos, endPos)
        while (i > 0 && endPos - 1 >= i && c[endPos - 1] != 0) // we do not check for i == 0 to avoid an underflow.
            --endPos;

        if (i != endPos)
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

            auto delegate = [&hits, &itExact, &itAll, bb, overlap, &params, &needles](
                TBiIter it, TNeedlesOverlap const & /*read*/, unsigned const errors_spent)
            {
                // TODO: we could turn reporting of exact iterators off at compile time by setting reportExactMatch = false if opt.directory is true. Evaluate binary size vs. performance.
                // WARNING: if it is computed on the directory, csvComputation currently still needs the exact matches (can be updated down below)
                if (errors_spent == 0)
                {
                    extend<true, csvComputation, errors>(it, hits, itExact, itAll, errors - errors_spent, needles, params.length,
                        params.length - overlap, params.length - 1, // searched interval
                        0, bb // entire interval
                    );
                }
                else
                {
                    extend<false, csvComputation, errors>(it, hits, itExact, itAll, errors - errors_spent, needles, params.length,
                        params.length - overlap, params.length - 1, // searched interval
                        0, bb // entire interval
                    );
                }
            };

            if (params.revCompl)
            {
                ModifiedString<ModifiedString<typename std::remove_reference<decltype(needles)>::type, ModComplementDna>, ModReverse> needlesRevCompl(needles);
                ModifiedString<ModifiedString<typename std::remove_reference<decltype(needlesOverlap)>::type, ModComplementDna>, ModReverse> needlesRevComplOverlap(needlesOverlap);
                using TNeedlesRevComplOverlap = decltype(needlesRevComplOverlap);

                // TODO: could store the exact hits as well and use these values!
                auto delegateRevCompl = [&hits, &itExact, &itAllrevCompl, bb, overlap, &params, &needlesRevCompl](
                    TBiIter it, TNeedlesRevComplOverlap const & /*read*/, unsigned const errors_spent)
                {
                    extend<false, csvComputation, errors>(it, hits, itExact, itAllrevCompl, errors - errors_spent, needlesRevCompl, params.length,
                        params.length - overlap, params.length - 1, // searched interval
                        0, bb // entire interval
                    );
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
                SEQAN_IF_CONSTEXPR (csvComputation)
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
                            distinct_sequences.insert(mappingSeqIdFile[location.i1]);
                        assert(entry.second.second.size() == 0 || params.revCompl);
                        for (auto const & location : entry.second.second) // reverse strand
                            distinct_sequences.insert(mappingSeqIdFile[location.i1]);

                        hits[j - beginPos] = distinct_sequences.size();

                        // NOTE: If you want to filter certain k-mers in the csv file based on the mappability value
                        // (with respect to --exclude-pseudo) you can unset 'entry' here.

                    }

                    if (!directory && countOccurrences(itExact[j - beginPos]) > 1)
                    {
                        // the for-loop does not insert an entry for kmers originating from a position such that the kmer spans two sequences. Hence we insert it here. The occurrences will later be cleared by resetLimits, but at least the position exists in the map.
                        myPosLocalize(entry.first, j, chromCumLengths); // TODO: inefficient for read data sets   0 > 0
                        if (entry.first.i2 > chromLengths[entry.first.i1] - params.length)
                        {
                            #pragma omp critical
                            locations.insert(entry);
                        }

                        for (auto const & exact_occ : getOccurrences(itExact[j - beginPos]))
                        {
                            entry.first = exact_occ;
                            // TODO: avoid copying
                            #pragma omp critical
                            locations.insert(entry);
                        }
                    }
                    else
                    {
                        myPosLocalize(entry.first, j, chromCumLengths); // TODO: inefficient for read data sets
                        // TODO: avoid copying
                        #pragma omp critical
                        locations.insert(entry);
                    }
                }

                if (!directory && countOccurrences(itExact[j - beginPos]) > 1) // guaranteed to exist, since there has to be at least one match!
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

        printProgress<outputProgress>(progressCount, progressStep, progressMax);
    }

    // The algorithm searches k-mers in the concatenation of all strings in the fasta file (e.g. chromosomes).
    // Hence, it also searches k-mers that overlap two strings that actually do not exist.
    // At the end we overwrite the frequency of those k-mers with 0.
    // TODO: k-mers spanning two strings should not be searched if there are many short strings (i.e., fasta of reads).
    resetLimits<csvComputation>(c, params.length, chromLengths, chromCumLengths, locations);
}
