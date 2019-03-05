#include <type_traits>

#include <seqan/arg_parse.h>

#include "genmap_helper.hpp"
#include "indexing.hpp"
#include "mappability.hpp"

template <typename TSpec, typename TLengthSum, unsigned LEVELS, unsigned WORDS_PER_BLOCK>
unsigned GemMapFastFMIndexConfig<TSpec, TLengthSum, LEVELS, WORDS_PER_BLOCK>::SAMPLING = 10;

using namespace seqan;

ArgumentParser::ParseResult parseCommandLineMain(int const argc, char const ** argv);

int main(int const argc, char const ** argv)
{
    if (std::string(CMAKE_BUILD_TYPE) != "Release")
    {
        std::cerr << "WARNING: This binary was not built in Release mode"
                     " and is expected to be much slower.\n";
    }

    // TODO: warning if no SSE4 and popcount support is used, that is might crash. use newer CPU

    // --version-check expects a parameter, the others (--copyright, --version) don't.
    int until = argc;
    bool skipNext = false;
    for (int i = 1; i < argc; ++i)
    {
        // version check expects a parameter
        if (std::string(argv[i]) == "--version-check")
            skipNext = true;

        if (argv[i][0] != '-')
        {
            if (skipNext)
            {
                skipNext = false;
            } else
            {
                until = i + 1;
                break;
            }
        }
    }

    ArgumentParser::ParseResult res = parseCommandLineMain(until, argv);

    if (res == ArgumentParser::PARSE_ERROR)
        return ArgumentParser::PARSE_ERROR;
    else if (res != ArgumentParser::PARSE_OK)
        return 0;

    --until; // undo the "+ 1" above

    if (std::string(argv[until]) == "map")
    {
        return mappabilityMain(argc - until, argv + until);
    }
    else if (std::string(argv[until]) == "index")
    {
        return indexMain(argc - until, argv + until);
    }
    else
    {
        // should not be reached
        std::cerr << "WRONG ARGUMENTS!\n";
        return -1;
    }

    return 0;
}

ArgumentParser::ParseResult parseCommandLineMain(int const argc, char const ** argv)
{
    ArgumentParser parser("GenMap");

    setShortDescription(parser, "Fast and Exact Computation of Genome Mappability");

    addUsageLine(parser, "[\\fIOPTIONS\\fP] COMMAND [\\fICOMMAND-OPTIONS\\fP]");
    sharedSetup(parser);

    addArgument(parser, ArgParseArgument(ArgParseArgument::STRING, "COMMAND"));
    setHelpText(parser, 0, "The sub-program to execute. See below.");
    setValidValues(parser, 0, "index map");

    addTextSection(parser, "Available commands");
    addText(parser, "\\fBindex  \\fP– Creates an index for mappability computation.");
    addText(parser, "\\fBmap  \\fP– Computes the mappability (requires a pre-built index).");
    addText(parser, "To view the help page for a specific command, simply run 'genmap command --help'.");

    return parse(parser, argc, argv);
}
