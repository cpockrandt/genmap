#include <type_traits>
#include <sys/types.h>
#include <dirent.h>

#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <seqan/index.h>

#include "common.hpp"

using namespace seqan;

struct IndexOptions
{
    CharString indexPath;
    uint64_t seqNumber;
    uint64_t maxSeqLength;
    uint64_t totalLength;
    unsigned sampling;
    bool directory;
    bool useRadix;
    bool verbose;
};

bool hasEnding(std::string const & fullString, std::string ending)
{
    if (fullString.length() >= ending.length())
        return (0 == fullString.compare (fullString.length() - ending.length(), ending.length(), ending));
    else
        return false;
}

template <typename TSeqNo, typename TSeqPos, typename TBWTLen,
          typename TString, typename TStringSetConfig, typename TRadixSortTag>
void buildIndex(StringSet<TString, TStringSetConfig> /*const*/ & chromosomes, IndexOptions const & options,
                TRadixSortTag const & /**/)
{
    using TText = StringSet<TString, Owner<ConcatDirect<SizeSpec_<TSeqNo, TSeqPos> > > > ;
    using TFMIndexConfig = TGemMapFastFMIndexConfig<TBWTLen>;
    TFMIndexConfig::SAMPLING = options.sampling;
    using TUniIndexConfig = seqan::FMIndex<TRadixSortTag, TFMIndexConfig>;

    TText chromosomesConcat(chromosomes);
    clear(chromosomes);

    {
        uint32_t const bwtDigits = std::numeric_limits<TBWTLen>::digits;
        uint32_t const seqNoDigits = std::numeric_limits<TSeqNo>::digits;
        uint32_t const seqPosDigits = std::numeric_limits<TSeqPos>::digits;

        // Print some size information on the index.
        if (options.verbose)
        {
            CharString const alphabet = std::is_same<typename Value<TString>::Type, Dna>::value ? "dna4" : "dna5";
            std::cout << "Index was constructed using " << alphabet << " alphabet.\n"
                         "The BWT is represented by " << bwtDigits << " bit values.\n"
                         "The sampled suffix array is represented by pairs of " << seqNoDigits <<
                         " and " << seqPosDigits << " bit values." << std::endl;
        }

        // Store index dimensions and alphabet type.
        StringSet<CharString, Owner<ConcatDirect<> > > info;
        uint32_t const alphabetSize = 4 + std::is_same<typename Value<TString>::Type, Dna5>::value;
        std::string const directoryFlag = options.directory ? "true" : "false";
        appendValue(info, "alphabet_size:" + std::to_string(alphabetSize));
        appendValue(info, "sa_dimensions_i1:" + std::to_string(seqNoDigits));
        appendValue(info, "sa_dimensions_i2:" + std::to_string(seqPosDigits));
        appendValue(info, "bwt_dimensions:" + std::to_string(bwtDigits));
        appendValue(info, "sampling_rate:" + std::to_string(options.sampling));
        appendValue(info, "fasta_directory:" + directoryFlag);
        save(info, toCString(std::string(toCString(options.indexPath)) + ".info"));
    }

    // TODO: do not use radix for small files
    // TODO: print helpful information if running out of memory.

    {
        Index<TText, TUniIndexConfig> fwdIndex(chromosomesConcat);
        if (std::is_same<TRadixSortTag, RadixSortSACreateTag>::value)
            indexCreateProgress(fwdIndex, FibreSALF());
        else
        {
            std::cout << "Create fwd Index ... " << std::flush;
            indexCreate(fwdIndex, FibreSALF());
            std::cout << "done!\n";
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
            std::cout << "Create fwd Index ... " << std::flush;
            indexCreate(fwdIndex, FibreSALF());
            std::cout << "done!\n";
        }
        clear(getFibre(getFibre(getFibre(fwdIndex, FibreSA()), FibreSparseString()), FibreValues()));
        clear(getFibre(getFibre(getFibre(fwdIndex, FibreSA()), FibreSparseString()), FibreIndicators()));
        saveRev(fwdIndex, toCString(std::string(toCString(options.indexPath)) + ".rev"));
    }
}

template <typename TString, typename TStringSetConfig, typename TRadixSortTag>
void buildIndex(StringSet<TString, TStringSetConfig> /*const*/ & chromosomes, IndexOptions const & options,
                TRadixSortTag const & /**/)
{
    constexpr uint64_t max16bitUnsignedValue = std::numeric_limits<uint16_t>::max();
    constexpr uint64_t max32bitUnsignedValue = std::numeric_limits<uint32_t>::max();

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

    addOption(parser, ArgParseOption("d", "directory", "Set this flag if --genome points to a directory of fasta files."));

    // TODO: describe both algorithms in terms of space consumption (disk and RAM)
    addOption(parser, ArgParseOption("A", "algorithm", "Algorithm for suffix array construction (needed for the FM index).", ArgParseArgument::INPUT_FILE, "IN"));
	setDefaultValue(parser, "algorithm", "radix");
	setValidValues(parser, "algorithm", "radix skew");

    addOption(parser, ArgParseOption("S", "sampling", "Sampling rate of suffix array", ArgParseArgument::INTEGER, "INT"));
	setDefaultValue(parser, "sampling", 10);
    setMaxValue(parser, "sampling", "64");
    setMinValue(parser, "sampling", "1");

    addOption(parser, ArgParseOption("v", "verbose", "Outputs some additional information on the constructed index."));

    addOption(parser, ArgParseOption("xA", "seqno", "Number of sequences", ArgParseArgument::INTEGER, "INT"));
    hideOption(parser, "seqno");

    addOption(parser, ArgParseOption("xB", "seqpos", "Max length of sequences", ArgParseArgument::INTEGER, "INT"));
    hideOption(parser, "seqpos");

    addOption(parser, ArgParseOption("xC", "bwtlen", "Total length of all sequences", ArgParseArgument::INTEGER, "INT"));
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
    getOptionValue(options.sampling, parser, "sampling");
    toLower(algorithm);
    options.directory = isSet(parser, "directory");
    options.useRadix = algorithm == "radix";
    options.verbose = isSet(parser, "verbose");

    // Check whether the output path exists and is writeable!
    if (fileExists(toCString(options.indexPath)))
    {
        std::cerr << "ERROR: The output directory for the index already exists at " << options.indexPath << '\n'
             << "Please remove it, or choose a different location.\n";
        return ArgumentParser::PARSE_ERROR;
    }
    else if (mkdir(toCString(options.indexPath), S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH))
    {
        std::cerr << "ERROR: Cannot create output directory at " << options.indexPath << '\n';
        return ArgumentParser::PARSE_ERROR;
    }

    // Append prefix name for indices.
    if (back(options.indexPath) != '/')
        options.indexPath += "/";
    options.indexPath += "index";

    // Read fasta input file(s)
    StringSet<Dna5String> chromosomes;
    if (options.directory)
    {
        uint32_t filesLoaded = 0;
        StringSet<CharString> directoryInformation;

        std::stringstream filenames;

        DIR * dirp = opendir(toCString(genomePath));
        struct dirent * dp;
        while ((dp = readdir(dirp)) != NULL)
        {
            std::string file(dp->d_name);
            if (hasEnding(file, ".fa") || hasEnding(file, ".fasta") || hasEnding(file, ".fastq"))
            {
                ++filesLoaded;
                filenames << file << '\n';
                std::string fullPath = toCString(genomePath);
                if (!hasEnding(fullPath, "/"))
                    fullPath += "/";
                fullPath += file;
                SeqFileIn seqFileIn(toCString(fullPath));

                StringSet<CharString, Owner<ConcatDirect<> > > ids;
                StringSet<Dna5String> chromosomes2;
                readRecords(ids, chromosomes2, seqFileIn);
                for (uint64_t i = 0; i < length(chromosomes2); ++i)
                {
                    std::string id = toCString(static_cast<CharString>(ids[i]));
                    std::string len = std::to_string(length(chromosomes2[i]));
                    // TODO: would need this for single fasta files as well!
                    // if (length(chromosomes2[i]) == 0)
                    // {
                    //     std::cerr << "WARNING: Found empty sequence in " << file << ". Sequence will be skipped.\n";
                    //     continue;
                    // }
                    appendValue(directoryInformation, file + ";" + len + ";" + id); // toCString(id.substr(0, id.find(" ")))
                    appendValue(chromosomes, chromosomes2[i]);
                }
                clear(ids);
            }
        }
        closedir(dirp);

        save(directoryInformation, toCString(std::string(toCString(options.indexPath)) + ".ids"));

        std::cout << filesLoaded << " fasta files have been loaded";
        if (options.verbose)
            std::cout << ":\n" << filenames.str();
        else
            std::cout << " (run with --verbose to list the files)" << '\n';
    }
    else
    {
        StringSet<CharString, Owner<ConcatDirect<> > > ids;
        SeqFileIn seqFileIn(toCString(genomePath));
        readRecords(ids, chromosomes, seqFileIn);
        if (options.verbose)
            std::cout << "Number of sequences in the fasta file: " << length(chromosomes) << '\n';

        // Store ids from fasta.
        save(ids, toCString(std::string(toCString(options.indexPath)) + ".ids"));
    }

    if (length(chromosomes) == 0)
    {
        // TODO: rmdir
        std::cerr << "ERROR: The fasta file seems to be empty.\n";
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
        options.seqNumber = (static_cast<uint64_t>(1) << seqno) - 2;
    }
    if (isSet(parser, "seqpos"))
    {
        uint64_t seqpos;
        getOptionValue(seqpos, parser, "seqpos");
        options.maxSeqLength = (static_cast<uint64_t>(1) << seqpos) - 2;
    }
    if (isSet(parser, "bwtlen"))
    {
        uint64_t bwtlen;
        getOptionValue(bwtlen, parser, "bwtlen");
        options.totalLength = (static_cast<uint64_t>(1) << bwtlen) - 2;
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

    std::cout << "Index created successfully.\n";

    return 0;
}
