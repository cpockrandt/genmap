#include <vector>
#include <string>
#include <cstdint>

// TODO: investigate performance of buffer sizes!
#define     BUFFER_SIZE     4*1024*1024 // 4 MB

using namespace seqan;

inline std::string get_output_path(Options const & opt, SearchParams const & searchParams)
{
    std::string path = std::string(toCString(opt.outputPath)) + "_" +
                       std::to_string(opt.errors) + "_" +
                       std::to_string(searchParams.length);
    return path;
}

template <typename T>
inline void saveRawFreq(std::vector<T> const & c, std::string const & output_path)
{
    std::ofstream outfile(output_path, std::ios::out | std::ios::binary);
    outfile.write((const char*) &c[0], c.size() * sizeof(T));
    outfile.close();
}

template <typename T>
inline void saveRawMap(std::vector<T> const & c, std::string const & output_path)
{
    char buffer[BUFFER_SIZE];
    std::ofstream outfile(output_path);
    outfile.rdbuf()->pubsetbuf(buffer, BUFFER_SIZE);

    for (T const & v : c)
    {
        float const f = 1.0 / static_cast<float>(v);
        outfile.write(reinterpret_cast<const char*>(&f), sizeof(float));
    }
    outfile.close();
}

template <typename T>
inline void saveTxtFreq(std::vector<T> const & c, std::string const & output_path)
{
    char buffer[BUFFER_SIZE];
    std::ofstream outfile(output_path, std::ios::out | std::ofstream::binary);
    outfile.rdbuf()->pubsetbuf(buffer, BUFFER_SIZE);

    copy(c.begin(), c.end(), (std::ostream_iterator<T>(outfile), std::ostream_iterator<int>(outfile, " ")));
    outfile.close();
}

template <typename T>
inline void saveTxtMap(std::vector<T> const & c, std::string const & output_path)
{
    char buffer[BUFFER_SIZE];
    std::ofstream outfile(output_path);
    outfile.rdbuf()->pubsetbuf(buffer, BUFFER_SIZE);

    for (T const & v : c)
    {
        float const f = 1.0 / static_cast<float>(v);
        outfile.write(reinterpret_cast<const char*>(&f), sizeof(float));
        // outfile << (1.0 / static_cast<float>(v)) << ' ';
    }
    outfile.close();
}
