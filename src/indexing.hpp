#include <type_traits>
#include <sys/types.h>
#include <dirent.h>
#include <cctype>

#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <seqan/index.h>

#include "common.hpp"
#include "seqan_libdivsufsort.h"

namespace seqan {
    // allow implicit conversion of non Dna5 characters to N instead of throwing an error
    template <>
    struct FastaIgnoreOrAssertFunctor_<Dna5> {
        typedef typename FastaIgnoreFunctor_<Dna5>::Type TIgnore;
        typedef TIgnore Type;
    };
}

using namespace seqan;

struct IndexOptions
{
    CharString indexPath;
    uint64_t seqNumber;
    uint64_t maxSeqLength;
    uint64_t totalLength;
    unsigned sampling;
    bool directory;
    bool useSkew;
    bool verbose;
};

void getFileNamesInDirectory(char * path, std::vector<std::pair<std::string, std::string>> & filenames, std::vector<std::string> const & fastaFileTypes)
{
    DIR * d = opendir(path);
    if (d == NULL)
        return;
    struct dirent * dir;
    while ((dir = readdir(d)) != NULL)
    {
        if (dir-> d_type != DT_DIR)
        {
            std::string const file(dir->d_name);
            std::string const fileExtension = file.substr(file.find_last_of('.') + 1);
            if (std::find(fastaFileTypes.begin(), fastaFileTypes.end(), fileExtension) != fastaFileTypes.end())
            {
                filenames.push_back({std::string(path) + "/", file});
            }
        }
        else if (dir -> d_type == DT_DIR && strcmp(dir->d_name, ".") != 0 && strcmp(dir->d_name, "..") != 0)
        {
            char d_path[4096];
            sprintf(d_path, "%s/%s", path, dir->d_name);
            getFileNamesInDirectory(d_path, filenames, fastaFileTypes);
        }
    }
    closedir(d);
}

inline std::string extractFileName(std::string const & path)
{
    auto const pos = path.find_last_of('/');
    if (pos == std::string::npos) // no slash found, i.e. file.fa
        return path;
    else // slash found, i.e., ./file.fa file.fa ../file.fa /path/to/file.fa
        return path.substr(pos + 1);
}

template <typename TAlgo, typename TSeqNo, typename TSeqPos, typename TBWTLen, typename TChromosomes>
void buildIndex(TChromosomes & chromosomes, IndexOptions const & options)
{
    using TString = typename Value<TChromosomes>::Type;
    using TAlphabet = typename Value<TString>::Type;
    using TText = StringSet<String<TAlphabet, Packed<> >, Owner<ConcatDirect<SizeSpec_<TSeqNo, TSeqPos> > > >; // here we tell it to pack it
    using TFMIndexConfig = TGemMapFastFMIndexConfig<TBWTLen>;
    using TUniIndexConfig = FMIndex<TAlgo, TFMIndexConfig>;
    TFMIndexConfig::SAMPLING = options.sampling;

    constexpr bool isDna5 = std::is_same<TAlphabet, Dna5>::value;

    TText chromosomesConcat(std::move(chromosomes)); // strings are getting packed
    clear(chromosomes); // reduce memory footprint

    {
        uint32_t const bwtDigits = std::numeric_limits<TBWTLen>::digits;
        uint32_t const seqNoDigits = std::numeric_limits<TSeqNo>::digits;
        uint32_t const seqPosDigits = std::numeric_limits<TSeqPos>::digits;

        // Print some size information on the index.
        if (options.verbose)
        {
            std::cout << "Index will be constructed using " << (isDna5 ? "dna5/rna5" : "dna4/rna4") << " alphabet.\n"
                         "- The BWT is represented by " << bwtDigits << " bit values.\n"
                         "- The sampled suffix array is represented by pairs of " << seqNoDigits <<
                         " and " << seqPosDigits << " bit values.\n";
        }

        // Store index dimensions and alphabet type.
        StringSet<CharString, Owner<ConcatDirect<> > > info;
        uint32_t const alphabetSize = 4 + isDna5;
        std::string const directoryFlag = options.directory ? "true" : "false";
        appendValue(info, "alphabet_size:" + std::to_string(alphabetSize));
        appendValue(info, "sa_dimensions_i1:" + std::to_string(seqNoDigits));
        appendValue(info, "sa_dimensions_i2:" + std::to_string(seqPosDigits));
        appendValue(info, "bwt_dimensions:" + std::to_string(bwtDigits));
        appendValue(info, "sampling_rate:" + std::to_string(options.sampling));
        appendValue(info, "fasta_directory:" + directoryFlag);
        appendValue(info, "packed_text:true");
        save(info, toCString(std::string(toCString(options.indexPath)) + ".info"));
    }

    {
        Index<TText, TUniIndexConfig> fwdIndex(chromosomesConcat);
        std::cout << "Create fwd Index ... " << std::flush;
        if (std::is_same<TAlgo, AlgoDivSufSortTag<int32_t> >::value || std::is_same<TAlgo, AlgoDivSufSortTag<int64_t> >::value)
        {
            indexCreate(fwdIndex, FibreSALF(), Fwd(), toCString(options.indexPath));
        }
        else
        {
            indexCreate(fwdIndex, FibreSALF());
            saveFwd(fwdIndex, toCString(options.indexPath));
        }
        std::cout << "done!\n";
    }

    reverse(chromosomesConcat, Serial());

    {
        Index<TText, TUniIndexConfig> bwdIndex(chromosomesConcat);
        std::cout << "Create bwd Index ... " << std::flush;
        if (std::is_same<TAlgo, AlgoDivSufSortTag<int32_t> >::value || std::is_same<TAlgo, AlgoDivSufSortTag<int64_t> >::value)
        {
            indexCreate(bwdIndex, FibreSALF(), Rev(), toCString(std::string(toCString(options.indexPath)) + ".rev"));
        }
        else
        {
            indexCreate(bwdIndex, FibreSALF());
            // TODO: clear not necessary?
            // clear(getFibre(getFibre(getFibre(bwdIndex, FibreSA()), FibreSparseString()), FibreValues()));
            // clear(getFibre(getFibre(getFibre(bwdIndex, FibreSA()), FibreSparseString()), FibreIndicators()));
            saveRev(bwdIndex, toCString(std::string(toCString(options.indexPath)) + ".rev"));
        }
        std::cout << "done!\n";
    }
}

template <typename TAlgo, typename TChromosomes>
void buildIndex(TChromosomes & chromosomes, IndexOptions const & options)
{
    constexpr uint64_t max16bitUnsignedValue = std::numeric_limits<uint16_t>::max();
    constexpr uint64_t max32bitUnsignedValue = std::numeric_limits<uint32_t>::max();

    // Analyze dimensions of the index needed.
    // NOTE: actually <= maxXXbitUnsignedValue+1 should be sufficient
    if (options.seqNumber <= max16bitUnsignedValue && options.maxSeqLength <= max32bitUnsignedValue)
    {
        if (options.totalLength <= max32bitUnsignedValue)
            buildIndex<TAlgo, uint16_t, uint32_t, uint32_t>(chromosomes, options); // e.g. human genome
        else
            buildIndex<TAlgo, uint16_t, uint32_t, uint64_t>(chromosomes, options); // e.g. barley genome
    }
    else if (options.seqNumber <= max32bitUnsignedValue && options.maxSeqLength <= max16bitUnsignedValue)
        buildIndex<TAlgo, uint32_t, uint16_t, uint64_t>(chromosomes, options); // e.g. read data set
    else
        buildIndex<TAlgo, uint64_t, uint64_t, uint64_t>(chromosomes, options); // anything else
}

template <typename TChromosomes>
int buildIndex(TChromosomes & chromosomes, IndexOptions const & options)
{
    if (options.useSkew)
    {
        std::cout << "Skew needs at least `25n` byes on disk. "
                     "Make sure that your temp directory has enough space/quota, or change it via `export TMPDIR=/somewhere/with/more/space`.\n"
                  << std::flush;
        buildIndex<Nothing>(chromosomes, options);
    }
    else
    {
        constexpr uint64_t max32bitSignedValue = std::numeric_limits<int32_t>::max();

        bool const divsufsort32bit = options.totalLength + options.seqNumber < max32bitSignedValue;

        if (divsufsort32bit)
            std::cout << "Input fits into int32_t (<2GB), algorithm will need about `6n` main memory.\n";
        else
            std::cout << "Input does not fit into int32_t (>2GB), algorithm will need about `10n` main memory.\n";

        std::cout << "It might be more or less depending on the number and length of the individual sequences.\n"
                     "If you are running out of memory, you can try to reduce the memory consumption a bit by inreasing `-S`, e.g., use `-S 20` (up to 64).\n"
                     "Although this will slow down the algorithm to compute the mappability.\n" << std::flush;

        if (divsufsort32bit)
            buildIndex<AlgoDivSufSortTag<int32_t> >(chromosomes, options);
        else
            buildIndex<AlgoDivSufSortTag<int64_t> >(chromosomes, options);
    }

    std::cout << "Index created successfully.\n";

    return 0;
}

template <typename TDirInfo, typename TChromosomes>
void readFasta(std::string const & fullPath, std::string const & file, TDirInfo & directoryInformation, TChromosomes & chromosomes)
{
    SeqFileIn seqFileIn(toCString(fullPath));
    StringSet<CharString, Owner<ConcatDirect<> > > ids, ids_short;
    std::vector<uint64_t> seq_len;

    typedef typename SeqFileBuffer_<TChromosomes, void>::Type TSeqBuffer;
    TSeqBuffer seqBuffer;

    // reuse the memory of context(seqFileIn).buffer for seqBuffer (which has a different type but same sizeof(Alphabet))
    swapPtr(seqBuffer.data_begin, context(seqFileIn).buffer[1].data_begin);
    swapPtr(seqBuffer.data_end, context(seqFileIn).buffer[1].data_end);
    seqBuffer.data_capacity = context(seqFileIn).buffer[1].data_capacity;

    uint64_t nbr_sequences = 0;
    while (!atEnd(seqFileIn))
    {
        readRecord(context(seqFileIn).buffer[0], seqBuffer, seqFileIn);
        // skip empty sequences
        if (length(seqBuffer) == 0)
        {
            continue;
        }
        ++nbr_sequences;
        appendValue(chromosomes, seqBuffer);
        seq_len.push_back(length(seqBuffer));

        CharString const & id = context(seqFileIn).buffer[0];
        appendValue(ids, id);
        // truncate id after first space
        uint32_t whitespace_pos = 0;
        while (whitespace_pos < length(id) && !std::isspace(static_cast<unsigned char>(id[whitespace_pos])))
        {
            ++whitespace_pos;
        }
        appendValue(ids_short, prefix(id, whitespace_pos));
    }

    swapPtr(seqBuffer.data_begin, context(seqFileIn).buffer[1].data_begin);
    swapPtr(seqBuffer.data_end, context(seqFileIn).buffer[1].data_end);
    context(seqFileIn).buffer[1].data_capacity = seqBuffer.data_capacity;
    seqBuffer.data_capacity = 0;

    if (nbr_sequences == 0)
    {
        std::cerr << "WARNING: The fasta file " << fullPath << " seems to be empty. Excluded from indexing.\n";
        return;
    }

    // if shortened ids are still unique, use them instead
    {
        StringSet<CharString> ids_short_copy(ids_short);
        std::sort(begin(ids_short_copy), end(ids_short_copy));
        if (std::unique(begin(ids_short_copy), end(ids_short_copy)) == end(ids_short_copy))
        {
            ids = ids_short;
        }
    }

    for (uint64_t i = 0; i < seq_len.size(); ++i)
    {
        std::string const id = toCString(static_cast<CharString>(ids[i]));
        std::string const len = std::to_string(seq_len[i]);

        appendValue(directoryInformation, file + ";" + len + ";" + id);
    }
}

int indexMain(int const argc, char const ** argv)
{
    // Argument Parser
    ArgumentParser parser("GenMap index");
    sharedSetup(parser);
    addDescription(parser, "Index creation. Only supports DNA and RNA (A, C, G, T/U, N). "
                           "Other characters will be converted to N.\n"
                           "Choose between the following index construction algorithms (-A / --algorithm):\n"
                           "* divsufsort (recommended, faster, needs about `6n` space in main memory/RAM, `10n` for sequences >2GB),\n"
                           "* skew (needs more than `25n` space on secondary memory/disk, i.e., in TMPDIR),\n"
                           "where `n` is the total number of bases in your fasta file(s).");

    // sorted in descending lexicographical order, since setValidValues() prints them in this order
    std::vector<std::string> const fastaFileTypes {"fsa", "fna", "fastq", "fasta", "fas", "faa", "fa"};
    std::string fastaFileTypesHelpString;
    for (uint8_t i = 0; i < fastaFileTypes.size() - 1; ++i)
        fastaFileTypesHelpString += '.' + fastaFileTypes[i] + ' ';
    fastaFileTypesHelpString += "and ." + fastaFileTypes.back();

    addOption(parser, ArgParseOption("F", "fasta-file", "Path to the fasta file.", ArgParseArgument::INPUT_FILE, "IN"));
    setValidValues(parser, "fasta-file", fastaFileTypes);

    addOption(parser, ArgParseOption("FD", "fasta-directory", "Path to the directory of fasta files "
        "(indexes all " + fastaFileTypesHelpString + " files in there, not including subdirectories).",
        ArgParseArgument::INPUT_FILE, "IN"));

    addOption(parser, ArgParseOption("I", "index", "Path to the index.", ArgParseArgument::OUTPUT_FILE, "OUT"));
    setRequired(parser, "index");

    addOption(parser, ArgParseOption("A", "algorithm", "Algorithm for suffix array construction "
        "(needed for the FM index).", ArgParseArgument::STRING, "TEXT"));
    setDefaultValue(parser, "algorithm", "divsufsort");
    setValidValues(parser, "algorithm", std::vector<std::string>{"divsufsort", "skew"});

    addOption(parser, ArgParseOption("S", "sampling", "Sampling rate of suffix array",
        ArgParseArgument::INTEGER, "INT"));
    setDefaultValue(parser, "sampling", 10);
    setMaxValue(parser, "sampling", "64");
    setMinValue(parser, "sampling", "1");

    addOption(parser, ArgParseOption("v", "verbose", "Outputs some additional information on the constructed index."));

    addOption(parser, ArgParseOption("xa", "seqno", "Number of sequences.", ArgParseArgument::INTEGER, "INT"));
    hideOption(parser, "seqno");

    addOption(parser, ArgParseOption("xb", "seqpos", "Max length of sequences.", ArgParseArgument::INTEGER, "INT"));
    hideOption(parser, "seqpos");

    addOption(parser, ArgParseOption("xc", "bwtlen", "Total length of all sequences.",
        ArgParseArgument::INTEGER, "INT"));
    hideOption(parser, "bwtlen");

    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    bool const isSetFastaFile = isSet(parser, "fasta-file");
    bool const isSetFastaDirectory = isSet(parser, "fasta-directory");

    if (isSetFastaFile && isSetFastaDirectory)
    {
        std::cerr << "ERROR: You can only use eiher --fasta-file or --fasta-directory, not both.\n";
        return ArgumentParser::PARSE_ERROR;
    }
    else if (!isSetFastaFile && !isSetFastaDirectory)
    {
        std::cerr << "ERROR: You forgot to specify --fasta-file or --fasta-directory.\n";
        return ArgumentParser::PARSE_ERROR;
    }

    // Retrieve input parameter
    IndexOptions options;
    CharString fastaPath, algorithm;
    getOptionValue(options.indexPath, parser, "index");
    getOptionValue(algorithm, parser, "algorithm");
    getOptionValue(options.sampling, parser, "sampling");
    toLower(algorithm);
    options.directory = isSetFastaDirectory;
    if (isSetFastaDirectory)
    {
        getOptionValue(fastaPath, parser, "fasta-directory");
        struct stat st;
        if (!(stat(toCString(fastaPath), &st) == 0 && S_ISDIR(st.st_mode)))
        {
            std::cerr << "ERROR: The fasta directory does not exist!\n";
            return ArgumentParser::PARSE_ERROR;
        }
    }
    else
    {
        getOptionValue(fastaPath, parser, "fasta-file");
        if (!fileExists(toCString(fastaPath)))
        {
            std::cerr << "ERROR: The fasta file does not exist!\n";
            return ArgumentParser::PARSE_ERROR;
        }
    }

    options.useSkew = algorithm == "skew";
    options.verbose = isSet(parser, "verbose");

    // Check whether the index path exists and is writeable!
    if (fileExists(toCString(options.indexPath)))
    {
        std::cerr << "ERROR: The directory for the index already exists at " << options.indexPath << '\n'
                  << "       Please remove it, or choose a different location.\n";
        return ArgumentParser::PARSE_ERROR;
    }
    else if (mkdir(toCString(options.indexPath), S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH))
    {
        std::cerr << "ERROR: Cannot create directory at " << options.indexPath << '\n';
        return ArgumentParser::PARSE_ERROR;
    }

    CharString const indexPathDir = options.indexPath;

    // Append prefix name for indices.
    if (back(options.indexPath) != '/')
        options.indexPath += '/';
    options.indexPath += "index";

    // Read fasta input file(s)
    StringSet<String<Dna5, Packed<> >, Owner<ConcatDirect<> > > chromosomes;
    StringSet<CharString, Owner<ConcatDirect<> > > directoryInformation;

    if (options.directory)
    {
        // find all fasta files in directory
        std::vector<std::pair<std::string, std::string> > filenames; // {path, filename}, e.g., {"/path/to/", "genome.fa"}
        getFileNamesInDirectory(toCString(fastaPath), filenames, fastaFileTypes);
        std::sort(filenames.begin(), filenames.end(), [](auto const & a, auto const & b) { return a.second < b.second; });

        // check for duplicate file names
        for (uint32_t i = 0 ; i < filenames.size() - 1; ++i)
        {
            if (filenames[i].second == filenames[i + 1].second)
            {
                rmdir(toCString(indexPathDir));
                std::cerr << "ERROR: At least two fasta files with the same filename found (this is not supported)! Please rename them and run again.\n";
                std::cerr << "       " << filenames[i].first << filenames[i].second << "!\n";
                std::cerr << "       " << filenames[i + 1].first << filenames[i + 1].second << "!\n";
                return ArgumentParser::PARSE_ERROR;
            }
        }

        for (auto const & file : filenames)
        {
            readFasta(file.first + file.second, file.second, directoryInformation, chromosomes);
        }

        if (length(chromosomes) == 0)
        {
            rmdir(toCString(indexPathDir));
            std::cerr << "ERROR: No (non-empty) fasta file found!\n";
            return ArgumentParser::PARSE_ERROR;
        }

        std::cout << filenames.size() << " fasta files have been loaded (run with --verbose to list the files):\n";
        if (options.verbose)
        {
            for (auto const & file : filenames)
            {
                std::cout << file.first << file.second << '\n';
            }
        }
    }
    else
    {
        std::string const file = extractFileName(toCString(fastaPath));
        readFasta(toCString(fastaPath), file, directoryInformation, chromosomes);
    }

    if (length(chromosomes) == 0)
    {
        rmdir(toCString(indexPathDir));
        std::cerr << "ERROR: There is no non-empty sequence in the fasta file(s).\n";
        return ArgumentParser::PARSE_ERROR;
    }

    save(directoryInformation, toCString(std::string(toCString(options.indexPath)) + ".ids"));

    // check whether it can be converted to Dna4 and analyze the data for determining the index dimensions later.
    bool canConvert = true;
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

    // overwrite index dimensions
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
        // Conversion to Dna4 alphabet since no Ns are in the sequences.
        // Unnecessary copy, replace with ModifiedString/View.
        // Not relevant for memory peaks though.
        StringSet<String<Dna, Packed<> >, Owner<ConcatDirect<> > > chromosomesDna4;
        move(chromosomesDna4, chromosomes);
        clear(chromosomes);
        return buildIndex(chromosomesDna4, options);
    }
    else
    {
        return buildIndex(chromosomes, options);
    }
}
