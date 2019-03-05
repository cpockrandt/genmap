#include <benchmark/benchmark.h>

#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <seqan/index.h>

static constexpr bool outputProgress = false;

// #include "../src/common.hpp"
// #include "../src/algo.hpp"

using namespace seqan;

// template <typename TSpec, typename TLengthSum, unsigned LEVELS, unsigned WORDS_PER_BLOCK>
// unsigned GemMapFastFMIndexConfig<TSpec, TLengthSum, LEVELS, WORDS_PER_BLOCK>::SAMPLING = 10;

void BM_TTEST(benchmark::State& state, uint64_t y)
{
    for (auto _ : state)
    {
        std::vector<float> xx;
        xx.resize(10000);
        // TODO
    }
}

BENCHMARK_CAPTURE(BM_TTEST, bench_name, 3)->Unit(benchmark::kMillisecond);

// BENCHMARK_MAIN();
int main(int argc, char** argv)
{
    // Argument parser
    // ArgumentParser parser("SearchSchemes - Benchmarking");
    // addDescription(parser,
    //     "App for creating the benchmark of Optimum Search Schemes from the paper.");
    //
    // addOption(parser, ArgParseOption("G", "genome", "Path to the indexed genome", ArgParseArgument::INPUT_FILE, "IN"));
	// setRequired(parser, "genome");
    //
    // addOption(parser, ArgParseOption("R", "reads", "Path to the reads", ArgParseArgument::INPUT_FILE, "IN"));
	// setValidValues(parser, "reads", "fa fasta fastq");
	// setRequired(parser, "reads");
    //
    // ArgumentParser::ParseResult res = parse(parser, argc, argv);
    // if (res != ArgumentParser::PARSE_OK)
    //     return res == ArgumentParser::PARSE_ERROR;
    //
    // // Retrieve input parameters
    // CharString indexPath, readsPath;
    // getOptionValue(indexPath, parser, "genome");
    // getOptionValue(readsPath, parser, "reads");
    //
    // open(fm_index, toCString(indexPath), OPEN_RDONLY);
    // StringSet<CharString> ids;
    // SeqFileIn seqFileIn(toCString(readsPath));
    // readRecords(ids, reads, seqFileIn);
    //
    // for (unsigned i = 1; i < length(reads); ++i)
    // {
    //     if (length(reads[i]) != length(reads[0]))
    //     {
    //         std::cerr << "ERROR: Not all reads have the same length." << std::endl;
    //         return 1;
    //     }
    // }

    ::benchmark::Initialize(&argc, argv);
    ::benchmark::RunSpecifiedBenchmarks();
}
