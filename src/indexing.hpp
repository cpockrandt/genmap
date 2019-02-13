#include <type_traits>

#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <seqan/index.h>

#include "common.hpp"

using namespace std;
using namespace seqan;

struct IndexOptions
{
    CharString indexPath;
    uint64_t seqNumber;
    uint64_t maxSeqLength;
    uint64_t totalLength;
    bool useRadix;
    bool verbose;
};

template <typename TSeqNo, typename TSeqPos, typename TBWTLen,
          typename TString, typename TStringSetConfig, typename TRadixSortTag>
void buildIndex(StringSet<TString, TStringSetConfig> /*const*/ & chromosomes, IndexOptions const & options,
                TRadixSortTag const & /**/)
{
    using TText = StringSet<TString, Owner<ConcatDirect<SizeSpec_<TSeqNo, TSeqPos> > > > ;
    using TFMIndexConfig = TGemMapFastFMIndexConfig<TBWTLen>;
    using TUniIndexConfig = seqan::FMIndex<TRadixSortTag, TFMIndexConfig>;

    TText chromosomesConcat(chromosomes);
    clear(chromosomes);

    {
        CharString alphabet = std::is_same<typename Value<TString>::Type, Dna>::value ? "dna4" : "dna5";
        uint32_t bwtDigits = std::numeric_limits<TBWTLen>::digits;
        uint32_t seqNoDigits = std::numeric_limits<TSeqNo>::digits;
        uint32_t seqPosDigits = std::numeric_limits<TSeqPos>::digits;

        // Print some size information on the index.
        if (options.verbose)
        {
            cout << "Index was constructed using " << alphabet << " alphabet.\n"
                    "The BWT is represented by " << bwtDigits << " bit values.\n"
                    "The sampled suffix array is represented by pairs of " << seqNoDigits << " and " << seqPosDigits <<
                    " bit values." << endl;
        }

        // Store index dimensions and alphabet type.
        CharString info;
        info += std::is_same<typename Value<TString>::Type, Dna>::value ? "dna4" : "dna5";
        info += "_" + std::to_string(seqNoDigits) + "_" + std::to_string(seqPosDigits) + "_" + std::to_string(bwtDigits);
        save(info, toCString(std::string(toCString(options.indexPath)) + ".info"));
    }

    // TODO: print helpful information if running out of memory.

    {
        Index<TText, TUniIndexConfig> fwdIndex(chromosomesConcat);
        if (std::is_same<TRadixSortTag, RadixSortSACreateTag>::value)
            indexCreateProgress(fwdIndex, FibreSALF());
        else
        {
            cout << "Create fwd Index ... " << flush;
            indexCreate(fwdIndex, FibreSALF());
            cout << "done!\n";
        }
        save(fwdIndex, toCString(options.indexPath));
    }

    {
        reverse(chromosomesConcat);
        Index<TText, TUniIndexConfig> fwdIndex(chromosomesConcat);
        if (std::is_same<TRadixSortTag, RadixSortSACreateTag>::value)
            indexCreateProgress(fwdIndex, FibreSALF());
        else
        {
            cout << "Create fwd Index ... " << flush;
            indexCreate(fwdIndex, FibreSALF());
            cout << "done!\n";
        }
        clear(getFibre(getFibre(getFibre(fwdIndex, FibreSA()), FibreSparseString()), FibreValues()));
        clear(getFibre(getFibre(getFibre(fwdIndex, FibreSA()), FibreSparseString()), FibreIndicators()));
        save(fwdIndex, toCString(std::string(toCString(options.indexPath)) + ".rev"));
    }
}

template <typename TString, typename TStringSetConfig, typename TRadixSortTag>
void buildIndex(StringSet<TString, TStringSetConfig> /*const*/ & chromosomes, IndexOptions const & options,
                TRadixSortTag const & /**/)
{
    constexpr uint64_t max16bitUnsignedValue = numeric_limits<uint16_t>::max();
    constexpr uint64_t max32bitUnsignedValue = numeric_limits<uint32_t>::max();

    // Analyze dimensions of the index needed.
    // NOTE: actually <= maxXXbitUnsignedValue+1 should be sufficient
    if (options.seqNumber <= max16bitUnsignedValue && options.maxSeqLength <= max32bitUnsignedValue)
    {
        if (options.totalLength <= max32bitUnsignedValue)
            buildIndex<uint16_t, uint32_t, uint32_t>(chromosomes, options, TRadixSortTag()); // e.g. human genome
        else
            buildIndex<uint16_t, uint32_t, uint64_t>(chromosomes, options, TRadixSortTag()); // e.g. barley genome
    }
    else if (options.seqNumber <= max32bitUnsignedValue && options.maxSeqLength <= max16bitUnsignedValue)
        buildIndex<uint32_t, uint16_t, uint64_t>(chromosomes, options, TRadixSortTag()); // e.g. read data set
    else
        buildIndex<uint64_t, uint64_t, uint64_t>(chromosomes, options, TRadixSortTag()); // anything else
}

template <typename TString, typename TStringSetConfig>
void buildIndex(StringSet<TString, TStringSetConfig> /*const*/ & chromosomes, IndexOptions const & options)
{
    if (options.useRadix)
        buildIndex(chromosomes, options, RadixSortSACreateTag());
    else
        buildIndex(chromosomes, options, Nothing());
}

int indexMain(int const argc, char const ** argv)
{
    // Argument Parser
    ArgumentParser parser("Index Creation");
    addDescription(parser, "App for creating an index. Only supports Dna (with and without N's).");

    addOption(parser, ArgParseOption("G", "genome", "Path to the genome", ArgParseArgument::INPUT_FILE, "IN"));
	setValidValues(parser, "genome", "fa fasta fastq");
	setRequired(parser, "genome");

    addOption(parser, ArgParseOption("I", "index", "Path to the index", ArgParseArgument::OUTPUT_FILE, "OUT"));
	setRequired(parser, "index");

    // TODO: describe both algorithms in terms of space consumption (disk and RAM)
    addOption(parser, ArgParseOption("A", "algorithm", "Algorithm for suffix array construction (needed for the FM index).", ArgParseArgument::INPUT_FILE, "IN"));
	setDefaultValue(parser, "algorithm", "radix");
	setValidValues(parser, "algorithm", "radix skew");

    addOption(parser, ArgParseOption("v", "verbose", "Outputs some additional information on the constructed index."));

    addOption(parser, ArgParseOption("a", "seqno", "Number of sequences", ArgParseArgument::INTEGER, "INT"));
    hideOption(parser, "seqno");

    addOption(parser, ArgParseOption("b", "seqpos", "Max length of sequences", ArgParseArgument::INTEGER, "INT"));
    hideOption(parser, "seqpos");

    addOption(parser, ArgParseOption("c", "bwtlen", "Total length of all sequences", ArgParseArgument::INTEGER, "INT"));
    hideOption(parser, "bwtlen");

    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    // Retrieve input parameter
    IndexOptions options;
    CharString genomePath, algorithm;
    getOptionValue(options.indexPath, parser, "index");
    getOptionValue(genomePath, parser, "genome");
    getOptionValue(algorithm, parser, "algorithm");
    toLower(algorithm);
    options.useRadix = algorithm == "radix";
    options.verbose = isSet(parser, "verbose");

    // Read fasta input file
    StringSet<CharString, Owner<ConcatDirect<> > > ids;
    StringSet<Dna5String> chromosomes;
    SeqFileIn seqFileIn(toCString(genomePath));
    readRecords(ids, chromosomes, seqFileIn);
    if (options.verbose)
        cout << "Number of sequences in the fasta file: " << length(chromosomes) << '\n';

    // Check whether the output path exists and is writeable!
    if (fileExists(toCString(options.indexPath)))
    {
        cerr << "ERROR: The output directory for the index already exists at " << options.indexPath << '\n'
             << "Please remove it, or choose a different location.\n";
        return ArgumentParser::PARSE_ERROR;
    }
    else if (mkdir(toCString(options.indexPath), S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH))
    {
        cerr << "ERROR: Cannot create output directory at " << options.indexPath << '\n';
        return ArgumentParser::PARSE_ERROR;
    }

    // Append prefix name for indices.
    if (back(options.indexPath) != '/')
        options.indexPath += "/";
    options.indexPath += "index";

    if (length(chromosomes) == 0)
    {
        cerr << "ERROR: The fasta file seems to be empty.\n";
        return ArgumentParser::PARSE_ERROR;
    }

    // check whether it can be converted to Dna4 and analyze the data for determining the index dimensions later.
    bool canConvert = true; // TODO: test this code block
    options.seqNumber = length(chromosomes);
    options.maxSeqLength = 0;
    options.totalLength = length(chromosomes); // to account for a sentinel character for each chromosome in the FM index.
    for (uint64_t i = 0; i < length(chromosomes); ++i)
    {
        options.totalLength += length(chromosomes[i]);
        options.maxSeqLength = std::max<uint64_t>(options.maxSeqLength, length(chromosomes[i]));

        for (uint64_t j = 0; canConvert && j < length(chromosomes[i]); ++j)
        {
            if (chromosomes[i][j] == 'N')
                canConvert = false;
        }
    }

    // overwrite dimensions
    if (isSet(parser, "seqno"))
    {
        uint64_t seqno;
        getOptionValue(seqno, parser, "seqno");
        options.seqNumber = (1ull << seqno) - 2;
    }
    if (isSet(parser, "seqpos"))
    {
        uint64_t seqpos;
        getOptionValue(seqpos, parser, "seqpos");
        options.maxSeqLength = (1ull << seqpos) - 2;
    }
    if (isSet(parser, "bwtlen"))
    {
        uint64_t bwtlen;
        getOptionValue(bwtlen, parser, "bwtlen");
        options.totalLength = (1ull << bwtlen) - 2;
    }

    // Construct index using Dna4 or Dna5 alphabet.
    if (canConvert)
    {
        // NOTE: avoid copying. Use a custom ModfiedFunctor instead.
        StringSet<DnaString> chromosomes4(chromosomes);
        clear(chromosomes);
        buildIndex(chromosomes4, options);
    }
    else
    {
        buildIndex(chromosomes, options);
    }

    // Store ids from fasta.
    save(ids, toCString(std::string(toCString(options.indexPath)) + ".ids"));

    cout << "Index created successfully.\n";

    return 0;
}
