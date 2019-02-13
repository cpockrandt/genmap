#include <chrono>

#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <seqan/index.h>

static constexpr bool outputProgress = false;

#include "common.h"
#include "algo.hpp"

using namespace std;
using namespace seqan;

template <typename TChar, typename TRng>
void randomText(String<TChar> & string, TRng & rng, uint64_t const length)
{
    resize(string, length);
    for (uint64_t i = 0; i < length; ++i)
        string[i] = TChar(rng() % 4);
}

template <unsigned errors, typename TIndex, typename TContainer>
inline void computeMappabilityTrivial(TIndex & index, TContainer & c, SearchParams const & searchParams)
{
    typedef typename TContainer::value_type value_type;

    constexpr uint64_t max_val = std::numeric_limits<typename TContainer::value_type>::max();
    auto scheme = OptimalSearchSchemes<0, errors>::VALUE;
    _optimalSearchSchemeComputeFixedBlocklength(scheme, searchParams.length);

    auto const & text = indexText(index);

    uint64_t global_pos = 0;
    for (uint64_t seq = 0; seq < length(text); ++seq)
    {
        for (uint64_t i = 0; i < length(text[seq]) - searchParams.length + 1; ++i, ++global_pos)
        {
            value_type hits = 0;
            auto delegate = [&hits](auto const &it, auto const & /*read*/, unsigned const /*errors*/) {
                if ((uint64_t) hits + countOccurrences(it) <= max_val)
                    hits += countOccurrences(it);
                else
                    hits = max_val;
            };

            auto const & needle = infix(text[seq], i, i + searchParams.length);
            Iter<TIndex, VSTree<TopDown<> > > it(index);
            _optimalSearchScheme(delegate, it, needle, scheme, HammingDistance());
            if (searchParams.revCompl)
            {
                DnaString needleRevCompl(needle);
                reverseComplement(needleRevCompl);
                goRoot(it);
                _optimalSearchScheme(delegate, it, needleRevCompl, scheme, HammingDistance());
            }
            c[global_pos] = hits;
        }
        global_pos += searchParams.length - 1;
    }
}

int main(int argc, char *argv[])
{
    using TIndexConfig = TBiIndexConfig<TGemMapFastFMIndexConfig<uint32_t, uint32_t>>;

    auto now = std::chrono::system_clock::now();
    auto seed = std::chrono::duration_cast<std::chrono::nanoseconds>(now.time_since_epoch()).count();
    cout << "Seed: " << seed << '\n';
    mt19937_64 rng(seed);

    constexpr unsigned errors = 3;

    while (true)
    {
        typedef StringSet<String<Dna>, Owner<ConcatDirect<> > > TGenome;
        TGenome genome;

        cout << "Stringset: ";
        uint64_t minStringLength = 999999999;
        uint64_t const stringsetsize = (rng() % 10) + 1;
        for (uint64_t ss = 0; ss < stringsetsize; ++ss)
        {
            uint64_t textLength = (rng() % 500) + 25;
            minStringLength = std::min(textLength, minStringLength);
            DnaString chr;
            randomText(chr, rng, textLength);
            appendValue(genome, chr);
            cout << textLength << ' ';
        }
        cout << endl;

        Index<TGenome, TIndexConfig> index(genome);
        indexCreate(index, FibreSALF());
        auto const & text = indexText(index).concat;

        for (uint64_t j = 0; j < 40; ++j)
        {
            uint64_t const length = std::min<uint64_t>((rand() % 25) + errors + 2, minStringLength);
            uint64_t const reads = (rand() % (length - 1 - errors)) + 1; // overlap must be at least errors+2
            uint64_t const overlap = length - reads + 1;

            cout << "Length (K): " << length << ", Overlap: " << overlap << ", Reads: " << reads << endl;
            vector<uint8_t> frequency_actual(seqan::length(text) - length + 1),
                            frequency_expected(seqan::length(text) - length + 1);

            SearchParams searchParams;
            searchParams.length = length;
            searchParams.overlap = overlap;
            searchParams.threads = omp_get_num_threads();
            searchParams.revCompl = rand() % 2;

            computeMappabilityTrivial<errors>(index, frequency_expected, searchParams);
            computeMappability<errors>(index, text, frequency_actual, searchParams);

            if (frequency_expected != frequency_actual)
            {
                cerr << "Length (K): " << length << ", Overlap: " << overlap << ", Reads: " << reads << endl;
                for (uint64_t ss = 0; ss < stringsetsize; ++ss)
                    cerr << genome[ss] << '\n';
                exit(1);
            }
        }
    }
}
