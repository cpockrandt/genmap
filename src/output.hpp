#include <vector>
#include <string>
#include <cstdint>

// TODO: investigate performance of buffer sizes!
// TODO: stack overflow might occur leading to a segmentation fault!
#define     BUFFER_SIZE     32*1024 // 32 KB

using namespace seqan;

inline std::string get_output_path(Options const & opt, SearchParams const & /*searchParams*/, std::string const & fastaFile)
{
    std::string path = std::string(toCString(opt.outputPath));
    if (back(path) != '/')
        path += "/";
    path += fastaFile;
                    // + fastaFile + "_" +
                    //    std::to_string(opt.errors) + "_" +
                    //    std::to_string(searchParams.length);
    return path;
}

// ---------------------------------------------------------------------------------------------------------------------

template <typename T>
void saveRawFreq(std::vector<T> const & c, std::string const & output_path)
{
    std::ofstream outfile(output_path, std::ios::out | std::ios::binary);
    outfile.write((const char*) &c[0], c.size() * sizeof(T));
    outfile.close();
}

template <typename T>
void saveRawMap(std::vector<T> const & c, std::string const & output_path)
{
    char buffer[BUFFER_SIZE];
    std::ofstream outfile(output_path);
    outfile.rdbuf()->pubsetbuf(buffer, BUFFER_SIZE);

    for (T const & v : c)
    {
        float const f = (v != 0) ? 1.0 / static_cast<float>(v) : 0;
        outfile.write(reinterpret_cast<const char*>(&f), sizeof(float));
    }
    outfile.close();
}

// ---------------------------------------------------------------------------------------------------------------------

template <typename T>
void saveTxtFreq(std::vector<T> const & c, std::string const & output_path)
{
    char buffer[BUFFER_SIZE];
    std::ofstream outfile(output_path + ".txt", std::ios::out | std::ofstream::binary);
    outfile.rdbuf()->pubsetbuf(buffer, BUFFER_SIZE);

    copy(c.begin(), c.end(), (std::ostream_iterator<T>(outfile), std::ostream_iterator<int>(outfile, " ")));
    outfile.close();
}

template <typename T>
void saveTxtMap(std::vector<T> const & c, std::string const & output_path)
{
    char buffer[BUFFER_SIZE];
    std::ofstream outfile(output_path + ".txt");
    outfile.rdbuf()->pubsetbuf(buffer, BUFFER_SIZE);

    for (T const & v : c)
    {
        float const f = (v != 0) ? 1.0 / static_cast<float>(v) : 0;
        // outfile.write(reinterpret_cast<const char*>(&f), sizeof(float));
        outfile << f << ' ';
    }
    outfile.close();
}

// ---------------------------------------------------------------------------------------------------------------------

template <bool mappability, typename T, typename TChromosomeNames, typename TChromosomeLengths>
void saveWig(std::vector<T> const & c, std::string const & output_path, TChromosomeNames const & chromNames, TChromosomeLengths const & chromLengths)
{
    uint64_t pos = 0;
    uint64_t begin_pos_string = 0;
    uint64_t end_pos_string = std::min(chromLengths[0], c.size());

    char buffer[BUFFER_SIZE];

    std::ofstream wigFile(output_path + ".wig");
    wigFile.rdbuf()->pubsetbuf(buffer, BUFFER_SIZE);

    for (uint64_t i = 0; i < length(chromLengths); ++i)
    {
        uint16_t current_val = c[pos];
        uint64_t occ = 0;
        uint64_t last_occ = 0;

        while (pos < end_pos_string)
        {
            if (current_val != c[pos])
            {
                if (last_occ != occ)
                    wigFile << "variableStep chrom=" << chromNames[i] << " span=" << occ << '\n';
                // TODO: document this behavior (mappability of 0)
                SEQAN_IF_CONSTEXPR (mappability)
                {
                    float const value = (current_val != 0) ? 1.0 / static_cast<float>(current_val) : 0;
                    wigFile << (pos - occ + 1 - begin_pos_string) << ' ' << value << '\n'; // pos in wig start at 1
                }
                else
                {
                    wigFile << (pos - occ + 1 - begin_pos_string) << ' ' << current_val << '\n'; // pos in wig start at 1
                }

                last_occ = occ;
                occ = 0;
                current_val = c[pos];
            }

            ++occ;
            ++pos;
        }

        // TODO: remove this block by appending a different value to c (reserve one more. check performance)
        if (last_occ != occ)
            wigFile << "variableStep chrom=" << chromNames[i] << " span=" << occ << '\n';
        SEQAN_IF_CONSTEXPR (mappability)
        {
            float const value = (current_val != 0) ? 1.0 / static_cast<float>(current_val) : 0;
            wigFile << (pos - occ + 1 - begin_pos_string) << ' ' << value << '\n'; // pos in wig start at 1
        }
        else
        {
            wigFile << (pos - occ + 1 - begin_pos_string) << ' ' << current_val << '\n'; // pos in wig start at 1
        }

        begin_pos_string += chromLengths[i];
        if (i + 1 < length(chromLengths))
            end_pos_string += chromLengths[i + 1];
        end_pos_string = std::min(end_pos_string, c.size()); // last chromosomeLength has to be reduced by K-1 characters
    }
    wigFile.close();

    // .chrom.sizes file
    std::ofstream chromSizesFile(output_path + ".chrom.sizes");
    for (uint64_t i = 0; i < length(chromLengths); ++i)
        chromSizesFile << chromNames[i] << '\t' << chromLengths[i] << '\n';
    chromSizesFile.close();

    // std::cout << "Wig file stored!" << std::endl;
}

// ---------------------------------------------------------------------------------------------------------------------

template <bool mappability, typename T, typename TChromosomeNames, typename TChromosomeLengths>
void saveBed(std::vector<T> const & c, std::string const & output_path, TChromosomeNames const & chromNames, TChromosomeLengths const & chromLengths)
{
    uint64_t pos = 0;
    uint64_t begin_pos_string = 0;
    uint64_t end_pos_string = std::min(chromLengths[0], c.size());

    char buffer[BUFFER_SIZE];

    std::ofstream bedFile(output_path + ".bed");
    bedFile.rdbuf()->pubsetbuf(buffer, BUFFER_SIZE);

    for (uint64_t i = 0; i < length(chromLengths); ++i)
    {
        uint16_t current_val = c[pos];
        uint64_t occ = 0;

        while (pos < end_pos_string)
        {
            if (current_val != c[pos])
            {
                bedFile << chromNames[i] << '\t'                    // chrom name
                        << (pos - occ - begin_pos_string) << '\t'   // start pos // pos in wig start at 1
                        << (pos - begin_pos_string - 1) << '\t'     // end pos
                        << '-' << '\t';                             // name

                SEQAN_IF_CONSTEXPR (mappability)
                    bedFile << ((current_val != 0) ? 1.0 / static_cast<float>(current_val) : 0) << '\n';
                else
                    bedFile << current_val << '\n';

                occ = 0;
                current_val = c[pos];
            }

            ++occ;
            ++pos;
        }

        // TODO: remove this block by appending a different value to c (reserve one more. check performance)
        bedFile << chromNames[i] << '\t'                    // chrom name
                << (pos - occ - begin_pos_string) << '\t'   // start pos // pos in wig start at 1
                << (pos - begin_pos_string - 1) << '\t'     // end pos
                << '-' << '\t';                             // name

        SEQAN_IF_CONSTEXPR (mappability)
            bedFile << ((current_val != 0) ? 1.0 / static_cast<float>(current_val) : 0) << '\n';
        else
            bedFile << current_val << '\n';

        begin_pos_string += chromLengths[i];
        if (i + 1 < length(chromLengths))
            end_pos_string += chromLengths[i + 1];
        end_pos_string = std::min(end_pos_string, c.size()); // last chromosomeLength has to be reduced by K-1 characters
    }
    bedFile.close();
}
