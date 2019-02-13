#include <vector>
#include <cstdint>
#include <limits>
#include <cmath>

#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <seqan/index.h>

using namespace std;
using namespace seqan;

static constexpr bool outputProgress = true;

#include "common.hpp"
#include "algo.hpp"

struct Options
{
    unsigned errors;
    bool mmap;
    bool indels;
    bool high;
    bool wigFile;
    CharString indexPath;
    CharString outputPath;
    CharString alphabet;
    uint32_t seqNoWidth;
    uint32_t maxSeqLengthWidth;
    uint32_t totalLengthWidth;
};

string get_output_path(Options const & opt, SearchParams const & searchParams)
{
    return string(toCString(opt.outputPath)) + "_" +
           to_string(opt.errors) + "_" +
           to_string(searchParams.length) + ".map" + (opt.high ? "16" : "8");
}

template <typename T>
inline void save(vector<T> const & c, string const & output_path)
{
    ofstream outfile(output_path, ios::out | ios::binary);
    outfile.write((const char*) &c[0], c.size() * sizeof(T));
    outfile.close();
}

template <typename TDistance, typename value_type, typename TIndex, typename TText>
inline void run(TIndex & index, TText const & text, Options const & opt, SearchParams const & searchParams)
{
    vector<value_type> c(length(text) - searchParams.length + 1, 0);

    switch (opt.errors)
    {
        case 0:  computeMappability<0>(index, text, c, searchParams);
                 break;
        case 1:  computeMappability<1>(index, text, c, searchParams);
                 break;
        case 2:  computeMappability<2>(index, text, c, searchParams);
                 break;
        case 3:  computeMappability<3>(index, text, c, searchParams);
                 break;
        case 4:  computeMappability<4>(index, text, c, searchParams);
                 break;
        default: cerr << "E = " << opt.errors << " not yet supported.\n";
                 exit(1);
    }

    if (outputProgress)
        std::cout << '\r';
    std::cout << "Progress: 100.00%\n" << std::flush;
    cout.flush();

    string output_path = get_output_path(opt, searchParams);
    save(c, output_path);

    if (opt.wigFile)
    {
        auto const & stringset = indexText(index);

        uint64_t pos = 0;
        uint64_t begin_pos_string = 0;
        uint64_t end_pos_string = length(stringset[0]);
        // for each sequence in the string set
        for (uint64_t i = 0; i < length(stringset); ++i)
        {
            stringstream ss;

            uint16_t current_val = c[pos];
            uint64_t occ = 0;
            uint64_t last_occ = 0;

            for (; pos < end_pos_string; ++pos)
            {
                if (current_val == c[pos])
                {
                    ++occ;
                }
                else
                {
                    if (last_occ != occ)
                        ss << "variableStep chrom=seq" << i << " span=" << occ << '\n';
                    float value = 1;
                    if (current_val != 0)
                        value = 1.0/(float)(current_val);
                    ss << (pos - occ + 1 - begin_pos_string) << ' ' << value << '\n'; // pos in wig start at 1

                    last_occ = occ;
                    occ = 1;
                    current_val = c[pos];
                }
            }

            if (last_occ != occ)
                ss << "variableStep chrom=seq" << i << " span=" << occ << '\n';
            float value = 1;
            if (current_val != 0)
                value = 1.0/(float)(current_val);
            ss << (pos - occ + 1 - begin_pos_string) << ' ' << value << '\n'; // pos in wig start at 1

            string wig_path = toCString(opt.outputPath);
            wig_path += "_" + to_string(opt.errors) + "_" + to_string(searchParams.length) + "_seq" + to_string(i);

            // .chrom.sizes file
            std::ofstream chromSizesFile;
            chromSizesFile.open(wig_path + ".chrom.sizes");
            chromSizesFile << "seq" << i << '\t' << length(stringset[i]) << '\n';
            chromSizesFile.close();

            // .wig file
            std::ofstream wigFile;
            wigFile.open(wig_path + ".wig");
            wigFile << ss.str();
            wigFile.close();

            begin_pos_string += length(stringset[i]);
            if (i + 1 < length(stringset))
                end_pos_string += length(stringset[i + 1]);
        }
    }
}

template <typename TChar, typename TAllocConfig, typename TDistance, typename value_type,
          typename TSeqNo, typename TSeqPos, typename TBWTLen>
inline void run(Options const & opt, SearchParams const & searchParams)
{
    typedef String<TChar, TAllocConfig> TString;
    typedef StringSet<TString, Owner<ConcatDirect<SizeSpec_<TSeqNo, TSeqPos> > > > TStringSet;

    using TFMIndexConfig = TGemMapFastFMIndexConfig<TBWTLen>;

    using TIndex = Index<TStringSet, TBiIndexConfig<TFMIndexConfig> >;

    TIndex index;
    open(index, toCString(opt.indexPath), OPEN_RDONLY);
    auto const & text = indexText(index);
    run<TDistance, value_type>(index, text.concat, opt, searchParams);
}

template <typename TChar, typename TAllocConfig, typename TDistance, typename TValue>
inline void run(Options const & opt, SearchParams const & searchParams)
{
    if (opt.seqNoWidth == 16 && opt.maxSeqLengthWidth == 32)
    {
        if (opt.totalLengthWidth == 32)
            run<TChar, TAllocConfig, TDistance, TValue, uint16_t, uint32_t, uint32_t>(opt, searchParams);
        else if (opt.totalLengthWidth == 64)
            run<TChar, TAllocConfig, TDistance, TValue, uint16_t, uint32_t, uint64_t>(opt, searchParams);
    }
    else if (opt.seqNoWidth == 32 && opt.maxSeqLengthWidth == 16 && opt.totalLengthWidth == 64)
        run<TChar, TAllocConfig, TDistance, TValue, uint32_t, uint16_t, uint64_t>(opt, searchParams);
    else if (opt.seqNoWidth == 64 && opt.maxSeqLengthWidth == 64 && opt.totalLengthWidth == 64)
        run<TChar, TAllocConfig, TDistance, TValue, uint64_t, uint64_t, uint64_t>(opt, searchParams);
}

template <typename TChar, typename TAllocConfig, typename TDistance>
inline void run(Options const & opt, SearchParams const & searchParams)
{
    if (opt.high)
        run<TChar, TAllocConfig, TDistance, uint16_t>(opt, searchParams);
    else
        run<TChar, TAllocConfig, TDistance, uint8_t>(opt, searchParams);
}

template <typename TChar, typename TAllocConfig>
inline void run(Options const & opt, SearchParams const & searchParams)
{
    if (opt.indels)
        run<TChar, TAllocConfig, EditDistance>(opt, searchParams);
    else
        run<TChar, TAllocConfig, HammingDistance>(opt, searchParams);
}

template <typename TChar>
inline void run(Options const & opt, SearchParams const & searchParams)
{
    if (opt.mmap)
        run<TChar, MMap<> >(opt, searchParams);
    else
        run<TChar, Alloc<> >(opt, searchParams);
}

int mappabilityMain(int argc, char const ** argv)
{
    // Argument parser
    ArgumentParser parser("GenMap");
    addDescription(parser,
        "Tool for computing the mappability on nucleotide sequences. It supports multi-fasta files with Dna4 and Dna5 alphabets.");

    addOption(parser, ArgParseOption("I", "index", "Path to the index", ArgParseArgument::INPUT_FILE, "IN"));
	setRequired(parser, "index");

    addOption(parser, ArgParseOption("O", "output", "Path to output directory (error number, length and overlap will be appended to the output file)", ArgParseArgument::OUTPUT_FILE, "OUT"));
    setRequired(parser, "output");

    addOption(parser, ArgParseOption("E", "errors", "Number of errors", ArgParseArgument::INTEGER, "INT"));

    addOption(parser, ArgParseOption("K", "length", "Length of k-mers", ArgParseArgument::INTEGER, "INT"));
    setRequired(parser, "length");

    addOption(parser, ArgParseOption("C", "reversecomplement", "Searches each k-mer on the reverse strand as well."));

    addOption(parser, ArgParseOption("i", "indels", "Turns on indels (EditDistance). "
        "If not selected, only mismatches will be considered."));

    addOption(parser, ArgParseOption("hi", "high", "Stores the mappability vector in 16 bit unsigned integers instead of 8 bit (max. value 65535 instead of 255)"));

    addOption(parser, ArgParseOption("w", "wig",
        "Output wig-files for adding a custom feature track to genome browsers. Mappability values will be stored as frequencies, i.e., 1/mappability. For each sequence (e.g., chromosome) a separate wig-file and chrom.size file is created"));

    addOption(parser, ArgParseOption("o", "overlap", "Number of overlapping reads (o + 1 Strings will be searched at once beginning with their overlap region). Default: K * (0.7^e * MIN(MAX(K,30),100) / 100)", ArgParseArgument::INTEGER, "INT"));
    //setRequired(parser, "overlap");

    addOption(parser, ArgParseOption("m", "mmap",
        "Turns memory-mapping on, i.e. the index is not loaded into RAM but accessed directly in secondary-memory. "
        "This makes the algorithm only slightly slower but the index does not have to be loaded into main memory "
        "(which takes some time)."));

    addOption(parser, ArgParseOption("t", "threads", "Number of threads", ArgParseArgument::INTEGER, "INT"));
    setDefaultValue(parser, "threads", omp_get_max_threads());

    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    // Retrieve input parameters
    Options opt;
    SearchParams searchParams;
    getOptionValue(opt.errors, parser, "errors");
    getOptionValue(opt.indexPath, parser, "index");
    getOptionValue(opt.outputPath, parser, "output");
    opt.mmap = isSet(parser, "mmap");
    opt.indels = isSet(parser, "indels");
    opt.high = isSet(parser, "high");
    opt.wigFile = isSet(parser, "wig");

    getOptionValue(searchParams.length, parser, "length");
    getOptionValue(searchParams.threads, parser, "threads");
    searchParams.revCompl = isSet(parser, "reversecomplement");

    if (opt.errors == 0)
        searchParams.overlap = searchParams.length * 0.7;
    else
        searchParams.overlap = searchParams.length * std::min(std::max(searchParams.length, 30u), 100u) * pow(0.7f, opt.errors) / 100.0;

    if (isSet(parser, "overlap"))
        getOptionValue(searchParams.overlap, parser, "overlap");
    // TODO: add verbose flag
    cout << "INFO: overlap = " << searchParams.overlap << '\n';

    if (searchParams.overlap > searchParams.length - 1)
    {
        cerr << "ERROR: overlap cannot be larger than K - 1.\n";
        return ArgumentParser::PARSE_ERROR;
    }

    if (!(searchParams.length - searchParams.overlap >= opt.errors + 2))
    {
        cerr << "ERROR: overlap should be at least K - E - 2. "
                "(K - O >= E + 2 must hold since common overlap has length K - O and will be split into E + 2 parts).\n";
        return ArgumentParser::PARSE_ERROR;
    }

    // searchParams.overlap - length of common overlap
    searchParams.overlap = searchParams.length - searchParams.overlap;

    if (opt.indels)
    {
        cerr << "ERROR: Indels are not supported yet.\n";
        return ArgumentParser::PARSE_ERROR;
    }

    // TODO: error message if output files already exist or directory is not writeable
    // TODO: nice error messages if index is incorrect or doesnt exist
    if (back(opt.indexPath) != '/')
        opt.indexPath += "/";
    opt.indexPath += "index";

    CharString info;
    open(info, toCString(std::string(toCString(opt.indexPath)) + ".info"));
    string infoStr(toCString(info));
    opt.alphabet = infoStr.substr(0, 4);
    opt.seqNoWidth = std::stoi(infoStr.substr(5, 2));
    opt.maxSeqLengthWidth = std::stoi(infoStr.substr(8, 2));
    opt.totalLengthWidth = std::stoi(infoStr.substr(11, 2));

    StringSet<CharString, Owner<ConcatDirect<> > > ids;
    open(ids, toCString(std::string(toCString(opt.indexPath)) + ".ids"));

    if (opt.alphabet == "dna4")
    {
        run<Dna>(opt, searchParams);
    }
    else
    {
        // run<Dna5>(opt, searchParams);
        std::cerr << "TODO: Dna5 alphabet has not been tested yet. Please do so and remove this error message afterwards.\n";
        exit(1);
    }

    return 0;
}
