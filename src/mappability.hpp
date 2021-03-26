#include <vector>
#include <cstdint>
#include <limits>
#include <sys/stat.h>
#include <set>

#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <seqan/bed_io.h>
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
    bool wigFile; // group files into mergable flags, i.e., BED | WIG, etc.
    bool bedFile;
    bool bedgraphFile;
    bool rawFile;
    bool txtFile;
    bool csvFile;
    bool designFile;
    bool outputPathIncludesFilename;
    OutputType outputType;
    bool directory;
    bool verbose;
    bool packed_text;
    CharString indexPath;
    CharString outputPath;
    CharString selectionPath;
    CharString alphabet;
    uint32_t seqNoWidth;
    uint32_t maxSeqLengthWidth;
    uint32_t totalLengthWidth;
    unsigned errors;
    unsigned sampling;
    uint64_t designWindowSize;
    // if rarest k-mer in window occurs < designPercentageSuperRare, pick all k-mers with this nbr. of hits
    // if rarest k-mer in window occurs < designPercentageRare (but not < SuperRare), pick only one k-mer
    float designPercentageSuperRare;
    float designPercentageRare;
};

struct DesignFileOutput
{
    // previously seen kmer -> kmer_id
    std::map<Dna5String, uint32_t> kmer_id;

    // file_id -> kmer_ids
    std::vector<std::set<uint32_t> > matrix;
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

    if (key == "packed_text") // this key was introduced later and might be missing in older indices
        return "false"; // older indices have unpacked/uncompressed texts

    // This should never happen unless the index file is corrupted or manipulated.
    std::cout << "ERROR: Malformed index.info file! Could not find key '" << key << "'.\n";
    exit(1);
}

template <typename TVector, typename TChromosomeNames, typename TChromosomeLengths, typename TLocations, typename TDirectoryInformation, typename TIntervals, typename TCSVIntervals, typename TText>
inline void outputMappability(TVector & c, Options const & opt, SearchParams const & searchParams,
                              std::string const & fastaFile, TChromosomeNames const & chromNames,
                              TChromosomeLengths const & chromLengths, TLocations & locations,
                              TDirectoryInformation const & directoryInformation,
                              TIntervals const & intervals, TCSVIntervals const & csvIntervals, bool const completeSameKmers,
                              DesignFileOutput & designFileOutput, uint64_t const currentFileNo, TText const & text,
                              TChromosomeLengths const & chromCumLengths)
{
    std::string output_path = std::string(toCString(opt.outputPath));
    if (!opt.outputPathIncludesFilename)
        output_path += fastaFile.substr(0, fastaFile.find_last_of('.')) + ".genmap";

    bool const outputSelection = opt.selectionPath != "";

    // reset mappability values that have been computed by accident using optimizations (copying values from same k-mers)
    if (outputSelection && completeSameKmers && (opt.rawFile || opt.txtFile || opt.wigFile || opt.bedgraphFile || opt.bedFile))
    {
        uint64_t last_interval_end = 0;
        for (auto const & interval : intervals) // triplets: chromosomeNamesId, interval.first, interval.second
        {
            for (uint64_t i = last_interval_end; i < std::get<0>(interval); ++i)
            {
                c[i] = 0;
            }
            last_interval_end = std::get<1>(interval);
        }

        for (uint64_t i = last_interval_end; i < c.size(); ++i)
        {
            c[i] = 0;
        }
    }

    if (opt.rawFile)
    {
        double start = get_wall_time();
        std::string output_path2 = output_path;
        if (opt.outputType == OutputType::mappability)
            output_path2 += ".map";
        else if (opt.outputType == OutputType::frequency_small)
            output_path2 += ".freq8";
        else // if (opt.outputType == OutputType::frequency_large)
            output_path2 += ".freq16";
        saveRaw(c, output_path2, (opt.outputType == OutputType::mappability));
        if (opt.verbose)
            std::cout << "- RAW file written in " << (round((get_wall_time() - start) * 100.0) / 100.0) << " seconds\n";
    }

    if (opt.txtFile)
    {
        double start = get_wall_time();
        saveTxt(c, output_path, chromNames, chromLengths, (opt.outputType == OutputType::mappability));
        if (opt.verbose)
            std::cout << "- TXT file written in " << (round((get_wall_time() - start) * 100.0) / 100.0) << " seconds\n";
    }

    if (opt.wigFile)
    {
        double start = get_wall_time();
        saveWig(c, output_path, chromNames, chromLengths, (opt.outputType == OutputType::mappability));
        if (opt.verbose)
            std::cout << "- WIG file written in " << (round((get_wall_time() - start) * 100.0) / 100.0) << " seconds\n";
    }

    if (opt.bedgraphFile)
    {
        double start = get_wall_time();
        saveBedGraph(c, output_path, chromNames, chromLengths, true /* bedgraph-file */, (opt.outputType == OutputType::mappability));
        if (opt.verbose)
            std::cout << "- bedgraph file written in " << (round((get_wall_time() - start) * 100.0) / 100.0) << " seconds\n";
    }

    if (opt.bedFile)
    {
        double start = get_wall_time();
        saveBedGraph(c, output_path, chromNames, chromLengths, false /* bed-file */, (opt.outputType == OutputType::mappability));
        if (opt.verbose)
            std::cout << "- BED file written in " << (round((get_wall_time() - start) * 100.0) / 100.0) << " seconds\n";
    }

    if (opt.csvFile)
    {
        double start = get_wall_time();
        saveCsv(output_path, locations, searchParams, directoryInformation, csvIntervals, outputSelection);
        if (opt.verbose)
            std::cout << "- CSV file written in " << (round((get_wall_time() - start) * 100.0) / 100.0) << " seconds\n";
    }

    if (opt.designFile)
    {
        double start = get_wall_time();
        if (opt.outputType == OutputType::mappability)
            saveDesignFile<true>(c, output_path, locations, searchParams, directoryInformation, csvIntervals, outputSelection, designFileOutput, currentFileNo, opt, text, chromCumLengths);
        else
            saveDesignFile<false>(c, output_path, locations, searchParams, directoryInformation, csvIntervals, outputSelection, designFileOutput, currentFileNo, opt, text, chromCumLengths);
        if (opt.verbose)
            std::cout << "- Design-File. file written in " << (round((get_wall_time() - start) * 100.0) / 100.0) << " seconds\n";
    }
}

template <typename TDistance, typename value_type, typename TSeqNo, typename TSeqPos,
          typename TIndex, typename TText, typename TChromosomeNames, typename TChromosomeLengths, typename TDirectoryInformation,
          typename TIntervals, typename TCSVIntervals>
inline void run(TIndex & index, TText const & text, Options const & opt, SearchParams const & searchParams,
                std::string const & fastaFile, TChromosomeNames const & chromNames, TChromosomeLengths const & chromLengths, TChromosomeLengths const & chromCumLengths,
                TDirectoryInformation const & directoryInformation, std::vector<TSeqNo> const & mappingSeqIdFile,
                TIntervals const & intervals, TCSVIntervals const & csvIntervals,
                uint64_t const currentFileNo, uint64_t const totalFileNo,
                DesignFileOutput & designFileOutput)
{
    std::vector<value_type> c(length(text), 0);

    std::map<Pair<TSeqNo, TSeqPos>,
             std::pair<std::vector<Pair<TSeqNo, TSeqPos> >,
                       std::vector<Pair<TSeqNo, TSeqPos> > > > locations;

    bool const csvComputation = opt.csvFile || searchParams.excludePseudo;
    bool completeSameKmers = true;

    switch (opt.errors)
    {
        case 0:  computeMappability<0>(index, text, c, searchParams, opt.directory, chromLengths, chromCumLengths, locations, mappingSeqIdFile, intervals, completeSameKmers, currentFileNo, totalFileNo, csvComputation);
                 break;
        case 1:  computeMappability<1>(index, text, c, searchParams, opt.directory, chromLengths, chromCumLengths, locations, mappingSeqIdFile, intervals, completeSameKmers, currentFileNo, totalFileNo, csvComputation);
                 break;
        case 2:  computeMappability<2>(index, text, c, searchParams, opt.directory, chromLengths, chromCumLengths, locations, mappingSeqIdFile, intervals, completeSameKmers, currentFileNo, totalFileNo, csvComputation);
                 break;
        case 3:  computeMappability<3>(index, text, c, searchParams, opt.directory, chromLengths, chromCumLengths, locations, mappingSeqIdFile, intervals, completeSameKmers, currentFileNo, totalFileNo, csvComputation);
                 break;
        case 4:  computeMappability<4>(index, text, c, searchParams, opt.directory, chromLengths, chromCumLengths, locations, mappingSeqIdFile, intervals, completeSameKmers, currentFileNo, totalFileNo, csvComputation);
                 break;
        default: std::cerr << "E > 4 not yet supported.\n";
                 exit(1);
    }
    SEQAN_IF_CONSTEXPR (outputProgress)
    {
        if (totalFileNo == 1)
        {
            std::cout << "\rProgress: 100.00%\x1b[K\n" << std::flush; // \e[K - clr_eol (remove anything after the cursor)
        }
        else
        {
            std::cout << "\r" // go up one line
                      << "File " << currentFileNo << " / " << totalFileNo << ". Progress: 100.00 %\x1b[K" << std::flush;

            if (opt.verbose || currentFileNo == totalFileNo)
            {
                std::cout << '\n'; // progress about writing files will follow
            }
        }
    }

    outputMappability(c, opt, searchParams, fastaFile, chromNames, chromLengths, locations, directoryInformation, intervals, csvIntervals, completeSameKmers, designFileOutput, currentFileNo, text, chromCumLengths);
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

    StringSet<CharString, Owner<ConcatDirect<> > > directoryInformation;
    open(directoryInformation, toCString(std::string(toCString(opt.indexPath)) + ".ids"), OPEN_RDONLY);
    appendValue(directoryInformation, "dummy.entry;0;chromosomename"); // dummy entry enforces that the mappability is
                                                                       // computed for the last file in the while loop.

    std::vector<TSeqNo> mappingSeqIdFile(length(directoryInformation) - 1);

    uint64_t totalFileNo;
    {
        uint64_t fastaId = 0;
        std::string fastaFile = std::get<0>(retrieveDirectoryInformationLine(directoryInformation[0]));
        for (uint64_t i = 0; i < length(directoryInformation) - 1; ++i)
        {
            auto const row = retrieveDirectoryInformationLine(directoryInformation[i]);
            if (std::get<0>(row) != fastaFile)
            {
                fastaFile = std::get<0>(row);
                ++fastaId;
            }
            if (searchParams.excludePseudo)
            {
                mappingSeqIdFile[i] = fastaId;
            }
        }
        totalFileNo = fastaId + 1;
    }

    // local begin and end positions (with respect to the corresponding sequence). id is sequence (std::string)
    std::map<std::string, std::vector<std::pair<uint64_t, uint64_t> > > intervals;
    if (opt.selectionPath != "")
    {
        BedFileIn bedIn(toCString(opt.selectionPath));
        BedRecord<Bed3> record;
        while (!atEnd(bedIn))
        {
            readRecord(record, bedIn);
            std::string const seqId{toCString(record.ref)};

            auto const lb = intervals.lower_bound(seqId);
            if(lb != intervals.end() && !(intervals.key_comp()(seqId, lb->first)))
                lb->second.emplace_back(std::make_pair(record.beginPos, record.endPos));
            else
                intervals.insert(lb, {seqId, {{record.beginPos, record.endPos}}});
        }
    }

    DesignFileOutput designFileOutput;
    designFileOutput.matrix.resize(totalFileNo);
    if (opt.designFile && opt.designPercentageRare < 1.0f / totalFileNo)
    {
        std::cerr << "There are only " << totalFileNo << " genomes, "
                  << "so the threshold Q=" << opt.designPercentageRare << " cannot be smaller than 1/" << totalFileNo << ".\n";
        exit(1);
    }

    if (opt.designFile && opt.designPercentageSuperRare < 1.0f / totalFileNo)
    {
        std::cerr << "WARNING: There are only " << totalFileNo << " genomes, "
                  << "so the threshold R=" << opt.designPercentageSuperRare << " cannot be smaller than 1/" << totalFileNo << " (unless you want to deactivate that feature).\n";
        // exit(1);
    }

    if (opt.designFile && opt.designPercentageRare < opt.designPercentageSuperRare)
    {
        std::cerr << "Threshold Q=" << opt.designPercentageRare << " should be greater or equal than R=" << opt.designPercentageSuperRare << "\n";
        exit(1);
    }

    auto const & text = indexText(index);
    std::map<std::string, uint64_t> chromosomeNamesDict;
    StringSet<CharString, Owner<ConcatDirect<> > > chromosomeNames;
    StringSet<uint64_t> chromosomeLengths, chromCumLengths; // ConcatDirect on PODs does not seem to support clear() ...
    uint64_t startPos = 0;
    uint64_t fastaFileLength = 0;
    uint64_t cumLength = 0;
    std::string fastaFile = std::get<0>(retrieveDirectoryInformationLine(directoryInformation[0]));
    appendValue(chromCumLengths, 0);
    uint64_t chromosomeNamesId = 0;
    // cumulative begin and end positions (with respect to the entire fasta file)
    std::vector<std::pair<uint64_t, uint64_t>> intervalsForSingleFasta;
    std::vector<std::tuple<uint32_t, uint64_t, uint64_t>> csvIntervalsForSingleFasta; // non-cumulative for csv-output filter

    std::vector<std::string> filenames;

    double start = get_wall_time();

    uint64_t currentFileNo = 0;

    std::vector<uint64_t> fasta_lengths;

    for (uint64_t i = 0; i < length(directoryInformation); ++i)
    {
        auto const row = retrieveDirectoryInformationLine(directoryInformation[i]);
        if (std::get<0>(row) != fastaFile)
        {
            //std::cout << "Now doing: " << fastaFile << std::endl;
            filenames.push_back(fastaFile);
            if (opt.csvFile || opt.designFile)
            {
                // sort for csv output
                std::sort(csvIntervalsForSingleFasta.begin(), csvIntervalsForSingleFasta.end(), [](auto const & t1, auto const & t2) {
                    if (std::get<0>(t1) != std::get<0>(t2))
                        return std::get<0>(t1) < std::get<0>(t2);
                    if (std::get<1>(t1) != std::get<1>(t2))
                        return std::get<1>(t1) < std::get<1>(t2);
                    return std::get<2>(t1) < std::get<2>(t2); // this is useless since the user should not input any overlapping intervals!
                });
            }

            ++currentFileNo;

            // do not compute/output mappability if only a subset shall be computed and the current fasta file does not contain any of the intervals of interest
            if (!(opt.selectionPath != "" && intervalsForSingleFasta.empty()))
            {
                // compute mappability for each fasta file
                auto const & fastaInfix = infixWithLength(text.concat, startPos, fastaFileLength);
                run<TDistance, value_type, TSeqNo, TSeqPos>(index, fastaInfix, opt, searchParams, fastaFile, chromosomeNames, chromosomeLengths, chromCumLengths, directoryInformation, mappingSeqIdFile, intervalsForSingleFasta, csvIntervalsForSingleFasta, currentFileNo, totalFileNo, designFileOutput);
            }

            fasta_lengths.push_back(fastaFileLength);
            startPos += fastaFileLength;
            fastaFile = std::get<0>(row);
            fastaFileLength = 0;
            chromosomeNamesDict.clear();
            clear(chromosomeNames);
            chromosomeNamesId = 0;
            clear(chromosomeLengths);
            clear(chromCumLengths);
            cumLength = 0;
            appendValue(chromCumLengths, 0);
            intervalsForSingleFasta.clear();
            csvIntervalsForSingleFasta.clear();
        }

        fastaFileLength += std::get<1>(row);
        chromosomeNamesDict.insert({toCString(std::get<2>(row)), chromosomeNamesId});

        // TODO: output warning if sequences in bed file do not appear in fasta/index
        auto intervalList = intervals.find(std::get<2>(row));
        if (intervalList != intervals.end())
        {
            for (auto const & interval : intervalList->second)
            {
                uint64_t const begin = cumLength + interval.first;
                uint64_t const end   = cumLength + interval.second;

                // check whether range is correct!
                if (interval.first >= std::get<1>(row) || interval.second > std::get<1>(row))
                {
                    std::cerr << "Error in BED file! Coordinates exceed sequence length: "
                              << "Seq. \"" << std::get<2>(row) << "\" has a length of " << std::get<1>(row) << ", "
                              << "but half-closed interval [" << interval.first << ", " << interval.second << ") given.\n";
                    exit(1);
                }

                intervalsForSingleFasta.emplace_back(std::make_pair(begin, end));

                if (opt.csvFile || opt.designFile)
                {
                    csvIntervalsForSingleFasta.emplace_back(std::make_tuple(chromosomeNamesId, interval.first, interval.second));
                }
            }
        }

        ++chromosomeNamesId;
        appendValue(chromosomeNames, std::get<2>(row));
        appendValue(chromosomeLengths, std::get<1>(row));
        cumLength += std::get<1>(row);
        appendValue(chromCumLengths, cumLength);
    }

    if (opt.designFile)
    {
        std::string output_path = std::string(toCString(opt.outputPath));

        uint32_t const nbr_total_kmers = designFileOutput.kmer_id.size();
        uint32_t probeCount = 0;
        uint32_t max_kmers_per_genome = 0;
        for (uint32_t i = 0; i < designFileOutput.matrix.size(); ++i)
        {
            // remove duplicates
            max_kmers_per_genome = std::max<uint64_t>(max_kmers_per_genome, designFileOutput.matrix[i].size());
            probeCount += designFileOutput.matrix[i].size();
        }

        // get median genome size
        std::sort(fasta_lengths.begin(), fasta_lengths.end());

        // output design file
        std::ofstream design_file(output_path + "genmap.nessie");
        design_file << "# Nessie database (strain-level classifier)\n";
        // TODO: add time and time zone
        design_file << "# build date: " << getDateTime() << '\n';
        design_file << "# window size: " << opt.designWindowSize << '\n';
        design_file << "# rare thresold: " << opt.designPercentageRare << '\n';
        design_file << "# super rare thresold: " << opt.designPercentageSuperRare << '\n';
        design_file << fasta_lengths[fasta_lengths.size()/2] << '\n';
        design_file << totalFileNo << '\t' << nbr_total_kmers << '\t' << max_kmers_per_genome << '\n';

        for (uint32_t i = 0; i < designFileOutput.matrix.size(); ++i)
        {
            design_file << filenames[i].substr(0, filenames[i].find_last_of(".")) << "\t1.0";

            for (const uint32_t k : designFileOutput.matrix[i])
            {
                design_file << '\t' << k;
            }
            for (uint32_t j = designFileOutput.matrix[i].size(); j < max_kmers_per_genome; ++j)
            {
                design_file << "\t0";
            }
            design_file << '\n';
        }

        // output kmer file (sorted by id)
        std::vector<Dna5String> kmer_vector;
        kmer_vector.resize(designFileOutput.kmer_id.size() + 1);
        for (auto & k : designFileOutput.kmer_id)
        {
            kmer_vector[k.second] = k.first;
        }

        design_file << "# kmers\n";

        for (uint32_t i = 1; i < kmer_vector.size(); ++i) // yes, it is correct that we start with 1 (see code block above)
        {
            design_file << i << '\t' << kmer_vector[i] << '\n';
        }

        design_file.close();
    }

    if (opt.verbose)
        std::cout << "Mappability computed in " << (round((get_wall_time() - start) * 100.0) / 100.0) << " seconds\n";
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

template <typename TChar>
inline void run(Options const & opt, SearchParams const & searchParams)
{
    if (opt.mmap && opt.packed_text)
        run<TChar, Packed<MMap<> >, HammingDistance>(opt, searchParams);
    else if (!opt.mmap && opt.packed_text)
        run<TChar, Packed<Alloc<> >, HammingDistance>(opt, searchParams);
    else if (opt.mmap && !opt.packed_text)
        run<TChar, MMap<>, HammingDistance>(opt, searchParams);
    else // if (!opt.mmap && !opt.packed_text)
        run<TChar, Alloc<>, HammingDistance>(opt, searchParams);
}

int mappabilityMain(int argc, char const ** argv)
{
    // Argument parser
    ArgumentParser parser("GenMap map");
    sharedSetup(parser);
    addDescription(parser,
        "Tool for computing the mappability/frequency on nucleotide sequences. It supports multi-fasta files with DNA or RNA alphabets (A, C, G, T/U, N). Frequency is the absolute number of occurrences, mappability is the inverse, i.e., 1 / frequency-value.");

    addOption(parser, ArgParseOption("I", "index", "Path to the index", ArgParseArgument::INPUT_FILE, "IN"));
	setRequired(parser, "index");

    addOption(parser, ArgParseOption("O", "output", "Path to output directory (or path to filename if only a single fasta files has been indexed)", ArgParseArgument::OUTPUT_FILE, "OUT"));
    setRequired(parser, "output");

    addOption(parser, ArgParseOption("E", "errors", "Number of errors", ArgParseArgument::INTEGER, "INT"));

    addOption(parser, ArgParseOption("K", "length", "Length of k-mers", ArgParseArgument::INTEGER, "INT"));
    setRequired(parser, "length");

    addOption(parser, ArgParseOption("S", "selection", "Path to a bed file (3 columns: chromosome, start, end) with selected coordinates to compute the mappability (e.g., exon coordinates)", ArgParseArgument::OUTPUT_FILE, "OUT"));

    addOption(parser, ArgParseOption("nc", "no-reverse-complement", "Searches the k-mers *NOT* on the reverse strand."));

    addOption(parser, ArgParseOption("ep", "exclude-pseudo", "Mappability only counts the number of fasta files that contain the k-mer, not the total number of occurrences (i.e., neglects so called- pseudo genes / sequences). This has no effect on the csv output."));

    addOption(parser, ArgParseOption("fs", "frequency-small", "Stores frequencies using 8 bit per value (max. value 255) instead of the mappbility using a float per value (32 bit). Applies to all formats (raw, txt, wig, bedgraph)."));
    addOption(parser, ArgParseOption("fl", "frequency-large", "Stores frequencies using 16 bit per value (max. value 65535) instead of the mappbility using a float per value (32 bit). Applies to all formats (raw, txt, wig, bedgraph)."));

    addOption(parser, ArgParseOption("r", "raw",
        "Output raw files, i.e., the binary format of std::vector<T> with T = float, uint8_t or uint16_t (depending on whether -fs or -fl is set). For each fasta file that was indexed a separate file is created. File type is .map, .freq8 or .freq16."));

    addOption(parser, ArgParseOption("t", "txt",
        "Output human readable text files, i.e., the mappability respectively frequency values separated by spaces (depending on whether -fs or -fl is set). For each fasta file that was indexed a separate txt file is created. WARNING: This output is significantly larger than raw files."));

    addOption(parser, ArgParseOption("w", "wig",
        "Output wig files, e.g., for adding a custom feature track to genome browsers. For each fasta file that was indexed a separate wig file and chrom.size file is created."));

    addOption(parser, ArgParseOption("bg", "bedgraph",
        "Output bedgraph files. For each fasta file that was indexed a separate bedgraph-file is created."));

    ArgParseOption bedFileOption("b", "bed",
        "Output bed files. For each fasta file that was indexed a separate bed-file is created.");
    hideOption(bedFileOption);
    addOption(parser, bedFileOption);

    addOption(parser, ArgParseOption("d", "csv",
                                     "Output a detailed csv file reporting the locations of each k-mer (WARNING: This will produce large files and makes computing the mappability significantly slower)."));

    addOption(parser, ArgParseOption("x", "design", "Output a design file for decoding strain abundance in a data set."));

    addOption(parser, ArgParseOption("W", "design-window", "Window size for k-mer extraction for design file", ArgParseArgument::INTEGER, "INT"));
    setDefaultValue(parser, "design-window", 1000);

    addOption(parser, ArgParseOption("Q", "design-rare-percentage", "Pick only one rare k-mer in window", ArgParseArgument::DOUBLE, "DOUBLE"));
    setDefaultValue(parser, "design-rare-percentage", 0.3);
    setMinValue(parser, "design-rare-percentage", "0.001");
    setMaxValue(parser, "design-rare-percentage", "1.0");

    addOption(parser, ArgParseOption("R", "design-super-rare-percentage", "Pick all super-rare k-mers in window", ArgParseArgument::DOUBLE, "DOUBLE"));
    setDefaultValue(parser, "design-super-rare-percentage", 0.1);
    setMinValue(parser, "design-super-rare-percentage", "0.001");
    setMaxValue(parser, "design-super-rare-percentage", "1.0");

    addOption(parser, ArgParseOption("m", "memory-mapping",
        "Turns memory-mapping on, i.e. the index is not loaded into RAM but accessed directly from secondary-memory. This may increase the overall running time, but do NOT use it if the index lies on network storage."));

    addOption(parser, ArgParseOption("T", "threads", "Number of threads", ArgParseArgument::INTEGER, "INT"));
    setDefaultValue(parser, "threads", omp_get_max_threads());

    addOption(parser, ArgParseOption("v", "verbose", "Outputs some additional information."));

    addOption(parser, ArgParseOption("xo", "overlap", "Number of overlapping reads (xo + 1 Strings will be searched at once beginning with their overlap region). Default: K * (0.7^e * MIN(MAX(K,30),100) / 100)", ArgParseArgument::INTEGER, "INT"));
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
    if (isSet(parser, "selection"))
        getOptionValue(opt.selectionPath, parser, "selection");

    opt.mmap = isSet(parser, "memory-mapping");
    opt.wigFile = isSet(parser, "wig");
    opt.bedgraphFile = isSet(parser, "bedgraph");
    opt.bedFile = isSet(parser, "bed");
    opt.rawFile = isSet(parser, "raw");
    opt.txtFile = isSet(parser, "txt");
    opt.csvFile = isSet(parser, "csv");
    opt.designFile = isSet(parser, "design");
    opt.verbose = isSet(parser, "verbose");

    getOptionValue(opt.designWindowSize, parser, "design-window");
    getOptionValue(opt.designPercentageRare, parser, "design-rare-percentage");
    getOptionValue(opt.designPercentageSuperRare, parser, "design-super-rare-percentage");

    if (!opt.wigFile && !opt.bedgraphFile && !opt.bedFile && !opt.rawFile && !opt.txtFile && !opt.csvFile && !opt.designFile)
    {
        std::cerr << "ERROR: Please choose at least one output format (i.e., --wig, --bedgraph, --bed, --raw, --txt, --csv).\n";
        return ArgumentParser::PARSE_ERROR;
    }

    // store in temporary variables to avoid parsing arguments twice
    bool const isSetFS = isSet(parser, "frequency-small");
    bool const isSetFL = isSet(parser, "frequency-large");

    if (isSetFS && isSetFL)
    {
        std::cerr << "ERROR: Cannot use both --frequency-small and --frequency-large. Please choose one.\n";
        return ArgumentParser::PARSE_ERROR;
    }

    if (isSetFS)
        opt.outputType = OutputType::frequency_small;
    else if (isSetFL)
        opt.outputType = OutputType::frequency_large;
    else // default value
        opt.outputType = OutputType::mappability;

    getOptionValue(searchParams.length, parser, "length");
    getOptionValue(searchParams.threads, parser, "threads");
    searchParams.revCompl = !isSet(parser, "no-reverse-complement");
    searchParams.excludePseudo = isSet(parser, "exclude-pseudo");

    // store in temporary variables to avoid parsing arguments twice
    bool const isSetOverlap = isSet(parser, "overlap");
    if (isSetOverlap)
        getOptionValue(searchParams.overlap, parser, "overlap");
    else if (opt.errors == 0)
        searchParams.overlap = searchParams.length * 0.7;
    else
        searchParams.overlap = searchParams.length * std::min(std::max(searchParams.length, 30u), 100u) * pow(0.7f, opt.errors) / 100.0;

    // (K - O >= E + 2 must hold since common overlap has length K - O and will be split into E + 2 parts)
    uint64_t const maxPossibleOverlap = std::min(searchParams.length - 1, searchParams.length - opt.errors - 2);
    if (searchParams.overlap > maxPossibleOverlap)
    {
        if (!isSetOverlap)
        {
            searchParams.overlap = maxPossibleOverlap;
        }
        else
        {
            std::cerr << "ERROR: overlap cannot be larger than min(K - 1, K - E - 2) = " << maxPossibleOverlap << ".\n";
            return ArgumentParser::PARSE_ERROR;
        }
    }

    // searchParams.overlap = length of common overlap
    searchParams.overlap = searchParams.length - searchParams.overlap;

    // TODO: error message if output files already exist or directory is not writeable
    // TODO: nice error messages if index is incorrect or doesnt exist
    if (back(opt.indexPath) != '/')
        opt.indexPath += '/';
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
    opt.packed_text = retrieve(info, "packed_text") == "true";

    // Check whether the output path exists
    {
        struct stat st;
        // is outputPath a directory and does it exist?
        if (stat(toCString(opt.outputPath), &st) == 0 && S_ISDIR(st.st_mode))
        {
            // okay (default case)
            opt.outputPathIncludesFilename = false;
            if (back(opt.outputPath) != '/')
                appendValue(opt.outputPath, '/');
        }
        // does outputPath include a filename?
        else if (!opt.directory)
        {
            // remove file name in temporary variable
            CharString outputPath2 = opt.outputPath;
            if (back(outputPath2) == '.')
            {
                // if it ends with . or .., it is always considered a directory
                appendValue(opt.outputPath, '/');
                opt.outputPathIncludesFilename = false;
            }
            else
            {
                int32_t last_slash_pos = length(outputPath2) - 1;
                // check for >= 0 in case it does not contain '/' at all (file in same directory)
                while (last_slash_pos >= 0 && outputPath2[last_slash_pos] != '/')
                    --last_slash_pos;
                if (last_slash_pos >= 0)
                {
                    erase(outputPath2, last_slash_pos, length(outputPath2));
                }
                else
                {
                    // since we checked at the very beginning whether it is an existing directory,
                    // we now assume that it is a filename and (in the current working directory)
                    // hence the path without the filename is '.'
                    outputPath2 = ".";
                }
                opt.outputPathIncludesFilename = true;
            }

            // check if the parent directory exists
            if (!(stat(toCString(outputPath2), &st) == 0 && S_ISDIR(st.st_mode)))
            {
                std::cerr << "ERROR: The output cannot be written to the file " << opt.outputPath << ".\n"
                          << "       It seems the directory " << outputPath2 << " does not exist.\n";
                return ArgumentParser::PARSE_ERROR;
            }
        }
        else
        {
            std::cerr << "ERROR: The output directory " << opt.outputPath << " does not exist.\n"
                      << "       A filename can only be specified for single indexed fasta files (not for indexed fasta directories).\n"
                      << "       Please create it, or choose a different location.\n";
            return ArgumentParser::PARSE_ERROR;
        }
    }

    if (opt.verbose)
    {
        // TODO: dna5/rna5
        std::cout << "Index was loaded (" << opt.alphabet << " alphabet, sampling rate of " << opt.sampling << ").\n"
                     "- The BWT is represented by " << opt.totalLengthWidth << " bit values.\n"
                     "- The sampled suffix array is represented by pairs of " << opt.seqNoWidth <<
                     " and " << opt.maxSeqLengthWidth  << " bit values.\n";

        if (opt.directory)
            std::cout << "- Index was built on an entire directory.\n" << std::flush;
        else
            std::cout << "- Index was built on a single fasta file.\n" << std::flush;
    }

    // TODO: remove opt.alphabet and replace by bool
    if (opt.alphabet == "dna4")
        run<Dna>(opt, searchParams);
    else
        run<Dna5>(opt, searchParams);

    return 0;
}
