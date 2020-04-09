#pragma once

#include <string>

#include <seqan/arg_parse.h>
#include <seqan/index.h>

using namespace seqan;

void sharedSetup(ArgumentParser & parser)
{
    // Set short description, version, and date.
    std::string versionString = SEQAN_APP_VERSION;
    setVersion(parser, versionString);
    setDate(parser, __DATE__);
    setShortCopyright(parser, "2019-2020 Christopher Pockrandt, released under the 3-clause-BSD; "
                              "2016-2019 Knut Reinert and Freie Universität Berlin, released under the 3-clause-BSD");

    setCitation(parser, "Pockrandt et al (2020). GenMap: Ultra-fast Computation of Genome Mappability.\n"
                        "doi: https://doi.org/10.1093/bioinformatics/btaa222");

    setLongCopyright(parser,
        " Copyright (c) 2019-2020, Christopher Pockrandt\n"
        " All rights reserved.\n"
        "\n"
        " This program is free software: you can redistribute it and/or modify\n"
        " it under the terms of the BSD-License (3-clause).\n"
        "\n"
        " GenMap is distributed in the hope that it will be useful,\n"
        " but WITHOUT ANY WARRANTY; without even the implied warranty of\n"
        " MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n"
        "\n"
        " You should have received a copy of the 3-clause BSD-License along with this\n"
        " program. If not, see <https://opensource.org/licenses/>.\n"
        "\n"
        " Copyright (c) 2013-2019 Hannes Hauswedell\n"
        " Copyright (c) 2016-2019 Knut Reinert and Freie Universität Berlin\n"
        " All rights reserved.\n"
        "\n"
        " Redistribution and use in source and binary forms, with or without\n"
        " modification, are permitted provided that the following conditions are met:\n"
        "\n"
        " * Redistributions of source code must retain the above copyright\n"
        "   notice, this list of conditions and the following disclaimer.\n"
        " * Redistributions in binary form must reproduce the above copyright\n"
        "   notice, this list of conditions and the following disclaimer in the\n"
        "   documentation and/or other materials provided with the distribution.\n"
        " * Neither the name of Knut Reinert or the FU Berlin nor the names of\n"
        "   its contributors may be used to endorse or promote products derived\n"
        "   from this software without specific prior written permission.\n"
        "\n"
        " THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS \"AS IS\"\n"
        " AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE\n"
        " IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE\n"
        " ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE\n"
        " FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL\n"
        " DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR\n"
        " SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER\n"
        " CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT\n"
        " LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY\n"
        " OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH\n"
        " DAMAGE.\n");

    addDescription(parser, "GenMap is a tool for fast and exact computation of genome mappability"
        " and can also be used for multiple genomes, e.g., to search for marker sequences.");

    addDescription(parser, "Detailed information is available in the wiki: "
        "<https://github.com/cpockrandt/genmap/wiki>");
}

template <typename TText, typename TSpec, typename TConfig>
inline bool open(Index<TText, BidirectionalIndex<FMIndex<TSpec, TConfig> > > & index, const char * fileName, int openMode)
{
    String<char> name;

    // fwd index
    name = fileName;    append(name, ".txt");
    if (!open(getFibre(index.fwd, FibreText()), toCString(name), openMode)) return false;

    name = fileName;    append(name, ".sa");
    if (!open(getFibre(index.fwd, FibreSA()), toCString(name), openMode)) return false;

    name = fileName;    append(name, ".lf");
    if (!open(getFibre(index.fwd, FibreLF()), toCString(name), openMode)) return false;

    setFibre(getFibre(index.fwd, FibreSA()), getFibre(index.fwd, FibreLF()), FibreLF());

    // rev index (only requires the BWT)
    name = fileName;    append(name, ".rev.lf");
    if (!open(getFibre(index.rev, FibreLF()), toCString(name), openMode)) return false;

    name = fileName;    append(name, ".sa.len");
    if (!open(index.rev.sa.sparseString._length, toCString(name), openMode)) return false;

    setFibre(getFibre(index.rev, FibreSA()), getFibre(index.rev, FibreLF()), FibreLF());

    return true;
}

template <typename TText, typename TSpec, typename TConfig>
inline bool saveRev(Index<TText, FMIndex<TSpec, TConfig> > const & index, const char * fileName, int openMode = OPEN_RDWR | OPEN_CREATE | OPEN_APPEND)
{
    String<char> name;

    // rev index (only requires the BWT)
    name = fileName;    append(name, ".lf");
    if (!save(getFibre(index, FibreLF()), toCString(name), openMode)) return false;

    return true;
}
