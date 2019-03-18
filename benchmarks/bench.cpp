#include <benchmark/benchmark.h>

#include <seqan/arg_parse.h>
#include <seqan/index.h>

static constexpr bool outputProgress = false;

#include "../src/genmap_helper.hpp"
#include "../src/mappability.hpp"

using namespace seqan;

// using TSeqNo = uint16_t;
// using TSeqPos = uint32_t;
// using TBWTLen = uint32_t;
//
// using TStringSet = StringSet<String<Dna, MMap<> >, Owner<ConcatDirect<SizeSpec_<TSeqNo, TSeqPos> > > >;
// using TFMIndexConfig = TGemMapFastFMIndexConfig<TBWTLen>;
// using TIndex = Index<TStringSet, TBiIndexConfig<TFMIndexConfig> >;
//
// using value_type = uint8_t;
//
// TIndex fmIndex;

template <typename TSpec, typename TLengthSum, unsigned LEVELS, unsigned WORDS_PER_BLOCK>
unsigned GemMapFastFMIndexConfig<TSpec, TLengthSum, LEVELS, WORDS_PER_BLOCK>::SAMPLING = 10;

seqan::CharString indexPath, outputPath;

template <uint64_t K, unsigned E>
void BM_GenMap(benchmark::State & state)
{
    constexpr bool csvComputation = false;

    Options opt;
    opt.mmap = true;
    opt.indels = false;
    opt.wigFile = false;
    opt.bedFile = false;
    opt.rawFile = true;
    opt.txtFile = false;
    opt.csvFile = false;
    opt.directory = false;
    opt.verbose = false;
    opt.outputType = OutputType::frequency_small;
    opt.indexPath = indexPath;
    opt.outputPath = outputPath;
    // CharString alphabet;
    // uint32_t seqNoWidth;
    // uint32_t maxSeqLengthWidth;
    // uint32_t totalLengthWidth;
    opt.errors = E;
    // unsigned sampling;

    SearchParams searchParams;
    searchParams.length = K;
    searchParams.overlap = computeOverlap(K, E, 0, false);
    searchParams.threads = 16;
    searchParams.revCompl = false;
    searchParams.excludePseudo = false;

    // std::vector<TSeqNo> mappingSeqIdFile; // only needed for excludePseudo = true
    // StringSet<uint64_t> chromosomeLengths;

    for (auto _ : state)
    {
        run(opt, searchParams);
        // auto const & text = indexText(fmIndex);
        // std::vector<value_type> c(lengthSum(text), 0);
        //
        // std::map<Pair<TSeqNo, TSeqPos>,
        //          std::pair<std::vector<Pair<TSeqNo, TSeqPos> >,
        //                    std::vector<Pair<TSeqNo, TSeqPos> > > > locations;
        //
        // computeMappability<E, csvComputation>(fmIndex, text.concat, c, searchParams, false /*directory*/, chromLengths, locations, mappingSeqIdFile);
    }
}

BENCHMARK_TEMPLATE(BM_GenMap, 36, 0)->Unit(benchmark::kMillisecond);
BENCHMARK_TEMPLATE(BM_GenMap, 24, 1)->Unit(benchmark::kMillisecond);
BENCHMARK_TEMPLATE(BM_GenMap, 36, 2)->Unit(benchmark::kMillisecond);
BENCHMARK_TEMPLATE(BM_GenMap, 50, 2)->Unit(benchmark::kMillisecond);
BENCHMARK_TEMPLATE(BM_GenMap, 75, 3)->Unit(benchmark::kMillisecond);

BENCHMARK_TEMPLATE(BM_GenMap, 101, 0)->Unit(benchmark::kMillisecond);
BENCHMARK_TEMPLATE(BM_GenMap, 101, 1)->Unit(benchmark::kMillisecond);
BENCHMARK_TEMPLATE(BM_GenMap, 101, 2)->Unit(benchmark::kMillisecond);
BENCHMARK_TEMPLATE(BM_GenMap, 101, 3)->Unit(benchmark::kMillisecond);
BENCHMARK_TEMPLATE(BM_GenMap, 101, 4)->Unit(benchmark::kMillisecond);

int main(int argc, char** argv)
{
    ArgumentParser parser("GenMap - Benchmarks");
    addDescription(parser, "Reproduces the benchmarks of GenMap published in the paper.");

    addOption(parser, ArgParseOption("I", "index", "Indexed genome", ArgParseArgument::INPUT_FILE, "IN"));
	setRequired(parser, "index");

    addOption(parser, ArgParseOption("O", "output", "Output directory", ArgParseArgument::OUTPUT_FILE, "OUT"));
	setRequired(parser, "output");

    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    // std::string indexPath;
    getOptionValue(indexPath, parser, "index");
    getOptionValue(outputPath, parser, "output");
    // if (back(indexPath) != '/')
    //     indexPath += '/';
    // indexPath += "index";

    // open(fmIndex, indexPath.c_str(), OPEN_RDONLY);
    // TFMIndexConfig::SAMPLING = 10;

    ::benchmark::Initialize(&argc, argv);
    ::benchmark::RunSpecifiedBenchmarks();
}
