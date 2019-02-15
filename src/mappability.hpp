#include <vector>
#include <cstdint>
#include <limits>
#include <cmath>

#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <seqan/index.h>

static constexpr bool outputProgress = true; // TODO: remove global variable

enum OutputType
{
    mappability,     // float (32 bit)
    frequency_large, // uint16_t (16 bit)
    frequency_small  // uint8_t (8 bit)
};

struct Options
{
    bool mmap;
    bool indels;
    bool wigFile; // group files into mergable flags, i.e., BED | WIG, etc.
    bool bedFile;
    bool rawFile;
    bool txtFile;
    OutputType outputType;
    bool directory;
    bool verbose;
    CharString indexPath;
    CharString outputPath;
    CharString alphabet;
    uint32_t seqNoWidth;
    uint32_t maxSeqLengthWidth;
    uint32_t totalLengthWidth;
    unsigned errors;
    unsigned sampling;
};

#include "common.hpp"
#include "algo.hpp"
#include "output.hpp"

using namespace seqan;

template <typename TSpec>
inline std::string retrieve(StringSet<CharString, TSpec> const & info, std::string const & key)
{
    for (uint32_t i = 0; i < length(info); ++i)
    {
        std::string row = toCString(static_cast<CharString>(info[i]));
        if (row.substr(0, length(key)) == key)
            return row.substr(length(key) + 1);
    }
    // This should never happen unless the index file is corrupted or manipulated.
    std::cout << "ERROR: Malformed index.info file! Could not find key '" << key << "'.\n";
    exit(1);
}

template <typename TDistance, typename value_type, typename TIndex, typename TText>
inline void run(TIndex & index, TText const & text, Options const & opt, SearchParams const & searchParams)
{
    std::vector<value_type> c(length(text) - searchParams.length + 1, 0);

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
        default: std::cerr << "E = " << opt.errors << " not yet supported.\n";
                 exit(1);
    }

    if (outputProgress)
        std::cout << '\r';
    std::cout << "Progress: 100.00%" << std::endl;

    std::string output_path = get_output_path(opt, searchParams);
    if (opt.rawFile)
    {
        if (opt.outputType == OutputType::mappability)
            saveRawMap(c, output_path + ".map");
        if (opt.outputType == OutputType::frequency_small)
            saveRawFreq(c, output_path + ".freq8");
        if (opt.outputType == OutputType::frequency_large)
            saveRawFreq(c, output_path + ".freq16");
    }

    if (opt.txtFile)
    {
        if (opt.outputType == OutputType::mappability)
            saveTxtMap(c, output_path + ".txt");
        if (opt.outputType == OutputType::frequency_small || opt.outputType == OutputType::frequency_large)
            saveTxtFreq(c, output_path + ".txt");
    }

    if (opt.bedFile)
    {
        std::cerr << "ERROR: bed file export is comming soon (in 1-2 days)!\n";
    }

    if (opt.wigFile)
    {
        auto const & stringset = indexText(index);

        uint64_t pos = 0;
        uint64_t begin_pos_string = 0;
        uint64_t end_pos_string = length(stringset[0]);
        // for each sequence in the string set
        for (uint64_t i = 0; i < length(stringset); ++i)
        {
            std::stringstream ss;

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

            std::string wig_path = toCString(opt.outputPath);
            wig_path += "_" + std::to_string(opt.errors) + "_" + std::to_string(searchParams.length) + "_seq" + std::to_string(i);

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
    TFMIndexConfig::SAMPLING = opt.sampling;

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
    if (opt.outputType == OutputType::frequency_large || opt.outputType == OutputType::mappability) // TODO: document precision for mappability
        run<TChar, TAllocConfig, TDistance, uint16_t>(opt, searchParams);
    else // if (opt.outputType == OutputType::frequency_small)
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
        "Tool for computing the mappability/frequency on nucleotide sequences. It supports multi-fasta files with Dna4 and Dna5 alphabets. Frequency is the absolute number of occurrences, mappability is the inverse, i.e., 1 / frequency-value.");

    addOption(parser, ArgParseOption("I", "index", "Path to the index", ArgParseArgument::INPUT_FILE, "IN"));
	setRequired(parser, "index");

    addOption(parser, ArgParseOption("O", "output", "Path to output directory and prefix filename (error number and length will be appended to the output file)", ArgParseArgument::OUTPUT_FILE, "OUT"));
    setRequired(parser, "output");

    addOption(parser, ArgParseOption("E", "errors", "Number of errors", ArgParseArgument::INTEGER, "INT"));

    addOption(parser, ArgParseOption("K", "length", "Length of k-mers", ArgParseArgument::INTEGER, "INT"));
    setRequired(parser, "length");

    addOption(parser, ArgParseOption("c", "reverse-complement", "Searches each k-mer on the reverse strand as well."));

    addOption(parser, ArgParseOption("i", "indels", "Turns on indels (EditDistance). "
        "If not selected, only mismatches will be considered."));

    addOption(parser, ArgParseOption("fs", "frequency-small", "Stores frequencies using 8 bit per value (max. value 255) instead of the mappbility using a float per value (32 bit). Applies to all formats (raw, txt, wig, bed)."));
    addOption(parser, ArgParseOption("fl", "frequency-large", "Stores frequencies using 16 bit per value (max. value 65535) instead of the mappbility using a float per value (32 bit). Applies to all formats (raw, txt, wig, bed)."));

    addOption(parser, ArgParseOption("r", "raw",
        "Output raw files, i.e., the binary format of std::vector<T> with T = float, uint8_t or uint16_t (depending on whether -fs or -fl is set). For each fasta file that was indexed a separate file is created. File type is .map, .freq8 or .freq16."));

    addOption(parser, ArgParseOption("t", "txt",
        "Output human readable text files, i.e., the mappability respectively frequency values separated by spaces (depending on whether -fs or -fl is set). For each fasta file that was indexed a separate txt file is created. WARNING: This output is significantly larger that raw files."));

    addOption(parser, ArgParseOption("w", "wig",
        "Output wig files, e.g., for adding a custom feature track to genome browsers. For each fasta file that was indexed a separate wig file and chrom.size file is created."));

    addOption(parser, ArgParseOption("b", "bed",
        "Output bed files. For each fasta file that was indexed a separate bed-file is created."));

    addOption(parser, ArgParseOption("m", "mmap",
        "Turns memory-mapping on, i.e. the index is not loaded into RAM but accessed directly in secondary-memory. "
        "This makes the algorithm only slightly slower but the index does not have to be loaded into main memory "
        "(which takes some time)."));

    addOption(parser, ArgParseOption("T", "threads", "Number of threads", ArgParseArgument::INTEGER, "INT"));
    setDefaultValue(parser, "threads", omp_get_max_threads());

    addOption(parser, ArgParseOption("v", "verbose", "Outputs some additional information."));

    addOption(parser, ArgParseOption("xo", "overlap", "Number of overlapping reads (o + 1 Strings will be searched at once beginning with their overlap region). Default: K * (0.7^e * MIN(MAX(K,30),100) / 100)", ArgParseArgument::INTEGER, "INT"));
    hideOption(parser, "overlap");

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
    opt.wigFile = isSet(parser, "wig");
    opt.bedFile = isSet(parser, "bed");
    opt.rawFile = isSet(parser, "raw");
    opt.txtFile = isSet(parser, "txt");
    opt.verbose = isSet(parser, "verbose");

    if (!opt.wigFile && !opt.bedFile && !opt.rawFile && !opt.txtFile)
    {
        std::cerr << "ERROR: Please choose at least one output format (i.e., --wig, --bed, --raw, --txt).\n";
        return ArgumentParser::PARSE_ERROR;
    }

    opt.outputType = OutputType::mappability; // default value
    if (isSet(parser, "frequency-small") && !isSet(parser, "frequency-large"))
    {
        opt.outputType = OutputType::frequency_small;
    }
    else if (!isSet(parser, "frequency-small") && isSet(parser, "frequency-large"))
    {
        opt.outputType = OutputType::frequency_large;
    }
    else if (isSet(parser, "frequency-small") && isSet(parser, "frequency-large"))
    {
        std::cerr << "ERROR: Cannot use both --frequency-small and --frequency-large. Please choose one.\n";
        return ArgumentParser::PARSE_ERROR;
    }


    getOptionValue(searchParams.length, parser, "length");
    getOptionValue(searchParams.threads, parser, "threads");
    searchParams.revCompl = isSet(parser, "reverse-complement");

    if (opt.errors == 0)
        searchParams.overlap = searchParams.length * 0.7;
    else
        searchParams.overlap = searchParams.length * std::min(std::max(searchParams.length, 30u), 100u) * pow(0.7f, opt.errors) / 100.0;

    if (isSet(parser, "overlap"))
        getOptionValue(searchParams.overlap, parser, "overlap");

    // std::cout << "INFO: overlap = " << searchParams.overlap << '\n';

    // TODO: make sure this never throws an error, if argument not set.
    if (searchParams.overlap > searchParams.length - 1)
    {
        std::cerr << "ERROR: overlap cannot be larger than K - 1.\n";
        return ArgumentParser::PARSE_ERROR;
    }

    if (!(searchParams.length - searchParams.overlap >= opt.errors + 2))
    {
        std::cerr << "ERROR: overlap should be at least K - E - 2. "
                "(K - O >= E + 2 must hold since common overlap has length K - O and will be split into E + 2 parts).\n";
        return ArgumentParser::PARSE_ERROR;
    }

    // searchParams.overlap - length of common overlap
    searchParams.overlap = searchParams.length - searchParams.overlap;

    // TODO: error message if output files already exist or directory is not writeable
    // TODO: nice error messages if index is incorrect or doesnt exist
    if (back(opt.indexPath) != '/')
        opt.indexPath += "/";
    opt.indexPath += "index";

    StringSet<CharString, Owner<ConcatDirect<> > > info;
    std::string infoPath = std::string(toCString(opt.indexPath)) + ".info";
    open(info, toCString(infoPath));
    opt.alphabet = "dna" + retrieve(info, "alphabet_size");
    opt.seqNoWidth = std::stoi(retrieve(info, "sa_dimensions_i1"));
    opt.maxSeqLengthWidth = std::stoi(retrieve(info, "sa_dimensions_i2"));
    opt.totalLengthWidth = std::stoi(retrieve(info, "bwt_dimensions"));
    opt.sampling = std::stoi(retrieve(info, "sampling_rate"));
    opt.directory = retrieve(info, "fasta_directory") == "true";

    StringSet<CharString, Owner<ConcatDirect<> > > ids;
    open(ids, toCString(std::string(toCString(opt.indexPath)) + ".ids"));

    if (opt.verbose)
    {
        std::cout << "Index was loaded (" << opt.alphabet << " alphabet, sampling rate of " << opt.sampling << ").\n"
                     "The BWT is represented by " << opt.totalLengthWidth << " bit values.\n"
                     "The sampled suffix array is represented by pairs of " << opt.seqNoWidth <<
                     " and " << opt.maxSeqLengthWidth  << " bit values.\n";

        if (opt.directory)
            std::cout << "Index was built on an entire directory (" << length(ids) << " fasta file(s))." << std::endl;
        else
            std::cout << "Index was built on a single fasta file." << std::endl;
    }

    if (opt.alphabet == "dna4")
    {
        run<Dna>(opt, searchParams);
    }
    else
    {
        std::cerr << "WARNING: Dna5 is still in beta-phase!\n";
        run<Dna5>(opt, searchParams);
    }

    return 0;
}
