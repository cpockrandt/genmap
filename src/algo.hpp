using namespace seqan;

template <typename TMappVector, typename TChromosomeLength>
void resetLimits(TMappVector const &, unsigned const, TChromosomeLength const)
{ }

template <typename TMappVector, typename TLengths, typename TConfig>
void resetLimits(TMappVector & c, unsigned const length, StringSet<TLengths, TConfig> const & chromLengths)
{
    for (unsigned i = 1; i < seqan::length(chromLengths) - 1; ++i)
    {
        for (unsigned j = 1; j < length; ++j)
        {
            c[chromLengths[i] - j] = 0;
        }
    }
}

// TODO: avoid signed integers

template <bool reportExactMatch, unsigned maxErrors, typename TBiIter, typename TValue, typename TText>
inline void extendExact(TBiIter it, TValue * hits, typename TBiIter::TFwdIndexIter * itExact,
                        TText const & text, unsigned const length,
                        uint64_t a, uint64_t b, // searched interval
                        uint64_t ab, uint64_t bb) // entire interval
{
    constexpr uint64_t max_val = std::numeric_limits<TValue>::max();

    if (b - a + 1 == length)
    {
        if (reportExactMatch && maxErrors == 0)
        {
            itExact[a-ab] = it.fwdIter;
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
                success = goDown(it2, text[i], Rev());
            }
            if (success)
                extendExact<reportExactMatch, maxErrors>(it2, hits, itExact, text, length, a, b_new, ab, bb);
        }
    //}

    if (a - 1 >= ab)
    {
        int64_t alm = b + 1 - length;
        uint64_t a_new = alm + std::max<int64_t>(((a - alm) - 1) >> 1, 0);
        for (int64_t i = a - 1; i >= static_cast<int64_t>(a_new); --i)
        {
            if(!goDown(it, text[i], Fwd()))
                return;
        }
        extendExact<reportExactMatch, maxErrors>(it, hits, itExact, text, length, a_new, b, ab, bb);
    }
}

// forward
template <bool reportExactMatch, unsigned maxErrors, typename TBiIter, typename TValue, typename TText>
inline void extend(TBiIter it, TValue * hits, typename TBiIter::TFwdIndexIter * itExact,
                   unsigned errorsLeft, TText const & text, unsigned const length,
                   uint64_t a, uint64_t b, // searched interval
                   uint64_t ab, uint64_t bb); // entire interval

template <bool reportExactMatch, unsigned maxErrors, typename TBiIter, typename TValue, typename TText>
inline void approxSearch(TBiIter it, TValue * hits, typename TBiIter::TFwdIndexIter * itExact,
                         unsigned errorsLeft, TText const & text, unsigned const length,
                         uint64_t a, uint64_t b, // searched interval
                         uint64_t ab, uint64_t bb, // entire interval
                         uint64_t b_new,
                         Rev const &)
{
    if (b == b_new)
    {
        extend<reportExactMatch, maxErrors>(it, hits, itExact, errorsLeft, text, length, a, b, ab, bb);
        return;
    }
    if (errorsLeft > 0)
    {
        if (goDown(it, Rev()))
        {
            do {
                bool delta = !ordEqual(parentEdgeLabel(it, Rev()), text[b + 1]);
                approxSearch<reportExactMatch, maxErrors>(it, hits, itExact, errorsLeft - delta, text, length, a, b + 1, ab, bb, b_new, Rev());
            } while (goRight(it, Rev()));
        }
    }
    else
    {
        for (uint64_t i = b + 1; i <= b_new; ++i)
        {
            if (!goDown(it, text[i], Rev()))
                return;
        }
        extendExact<reportExactMatch, maxErrors>(it, hits, itExact, text, length, a, b_new, ab, bb);
    }
}
template <bool reportExactMatch, unsigned maxErrors, typename TBiIter, typename TValue, typename TText>
inline void approxSearch(TBiIter it, TValue * hits, typename TBiIter::TFwdIndexIter * itExact,
                         unsigned errorsLeft, TText const & text, unsigned const length,
                         uint64_t a, uint64_t b, // searched interval
                         uint64_t ab, uint64_t bb, // entire interval
                         uint64_t a_new,
                         Fwd const & /*tag*/)
{
    if (a == a_new)
    {
        extend<reportExactMatch, maxErrors>(it, hits, itExact, errorsLeft, text, length, a, b, ab, bb);
        return;
    }
    if (errorsLeft > 0)
    {
        if (goDown(it, Fwd()))
        {
            do {
                bool delta = !ordEqual(parentEdgeLabel(it, Fwd()), text[a - 1]);
                approxSearch<reportExactMatch, maxErrors>(it, hits, itExact, errorsLeft - delta, text, length, a - 1, b, ab, bb, a_new, Fwd());
            } while (goRight(it, Fwd()));
        }
    }
    else
    {
        for (int64_t i = a - 1; i >= static_cast<int64_t>(a_new); --i)
        {
            if (!goDown(it, text[i], Fwd()))
                return;
        }
        extendExact<reportExactMatch, maxErrors>(it, hits, itExact, text, length, a_new, b, ab, bb);
    }
}

template <bool reportExactMatch, unsigned maxErrors, typename TBiIter, typename TValue, typename TText>
inline void extend(TBiIter it, TValue * hits, typename TBiIter::TFwdIndexIter * itExact,
                   unsigned errorsLeft, TText const & text, unsigned const length,
                   uint64_t a, uint64_t b, // searched interval
                   uint64_t ab, uint64_t bb) // entire interval
{
    constexpr uint64_t max_val = std::numeric_limits<TValue>::max();

    if (errorsLeft == 0)
    {
        extendExact<reportExactMatch, maxErrors>(it, hits, itExact, text, length, a, b, ab, bb);
        return;
    }
    if (b - a + 1 == length)
    {
        if (reportExactMatch && maxErrors == errorsLeft)
            itExact[a-ab] = it.fwdIter;
        hits[a-ab] = std::min((uint64_t) countOccurrences(it) + hits[a-ab], max_val);
        return;
    }
    //if (b + 1 <= bb)
    //{
        uint64_t brm = a + length - 1;
        uint64_t b_new = b + (((brm - b) + 2 - 1) >> 1); // ceil((bb - b)/2)
        if (b_new <= bb)
        {
            approxSearch<reportExactMatch, maxErrors>(it, hits, itExact, errorsLeft, text, length,
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
        approxSearch<reportExactMatch, maxErrors>(it, hits, itExact, errorsLeft, text, length,
                     a, b, // searched interval
                     ab, bb, // entire interval
                     a_new,
                     Fwd()
        );
    }
}

template <unsigned errors, typename TIndex, typename TText, typename TContainer, typename TChromosomeLengths>
inline void computeMappability(TIndex & index, TText const & text, TContainer & c, SearchParams const & params, bool const directory, TChromosomeLengths const & chromLengths)
{
    typedef typename TContainer::value_type TValue;
    typedef Iter<TIndex, VSTree<TopDown<> > > TBiIter;

    auto const & limits = stringSetLimits(indexText(index));
    uint64_t const textLength = length(text);
    uint64_t const numberOfKmers = textLength - params.length + 1;
    uint64_t const stepSize = params.length - params.overlap + 1; // Number of overlapping k-mers searched at once

    // Number of loop iterations assigned to a thread at once
    // It should be significantly smaller than the number of loop iterations (numberOfKmers/stepSize), since
    // the running time of different loop iterations can vary vastly (e.g., repeats are slower than unique regions).
    uint64_t const chunkSize = std::max<uint64_t>(1, numberOfKmers / (stepSize * params.threads * 50));

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

            auto scheme = OptimalSearchSchemes<0, errors>::VALUE;
            _optimalSearchSchemeComputeFixedBlocklength(scheme, overlap);

            typename TBiIter::TFwdIndexIter itExact[endPos - beginPos];
            TValue hits[endPos - beginPos], hitsRevCompl[endPos - beginPos];

            auto const & needles = infix(text, beginPos, beginPos + params.length + (endPos - beginPos) - 1);
            auto const & needlesOverlap = infix(text, beginPos + params.length - overlap, beginPos + params.length);
            using TNeedlesOverlap = decltype(needlesOverlap);

            uint64_t const bb = std::min(textLength - 1, params.length - 1 + params.length - overlap);

            auto delegate = [&hits, &itExact, bb, overlap, &params, &needles](
                TBiIter it, TNeedlesOverlap const & /*read*/, unsigned const errors_spent)
            {
                if (errors_spent == 0)
                {
                    extend<true, errors>(it, hits, itExact, errors - errors_spent, needles, params.length,
                        params.length - overlap, params.length - 1, // searched interval
                        0, bb // entire interval
                    );
                }
                else
                {
                    extend<false, errors>(it, hits, itExact, errors - errors_spent, needles, params.length,
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

                // uint64_t const bb = std::min(textLength - 1, params.length - 1 + params.length - overlap);

                auto delegateRevCompl = [&hitsRevCompl, &itExact, bb, overlap, &params, &needlesRevCompl](
                    TBiIter it, TNeedlesRevComplOverlap const & /*read*/, unsigned const errors_spent)
                {
                    extend<false, errors>(it, hitsRevCompl, itExact, errors - errors_spent, needlesRevCompl, params.length,
                        params.length - overlap, params.length - 1, // searched interval
                        0, bb // entire interval
                    );
                };

                for (unsigned y = 0; y < (endPos - beginPos); ++y)
                    hitsRevCompl[y] = 0;

                TBiIter it(index);
                _optimalSearchScheme(delegateRevCompl, it, needlesRevComplOverlap, scheme, HammingDistance());

                for (unsigned y = 0; y < (endPos - beginPos); ++y)
                    hits[(endPos - beginPos - 1) - y] = hitsRevCompl[y];
            }
            else
            {
                for (unsigned y = 0; y < (endPos - beginPos); ++y)
                    hits[y] = 0;
            }

            TBiIter it(index);
            _optimalSearchScheme(delegate, it, needlesOverlap, scheme, HammingDistance());
            for (uint64_t j = beginPos; j < endPos; ++j)
            {
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
    resetLimits(c, params.length, chromLengths);
}
