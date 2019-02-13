// ==========================================================================
//                                  lambda
// ==========================================================================
// Copyright (c) 2013-2017, Hannes Hauswedell <h2 @ fsfe.org>
// Copyright (c) 2016-2017, Knut Reinert and Freie Universität Berlin
// All rights reserved.
//
// This file is part of Lambda.
//
// Lambda is Free Software: you can redistribute it and/or modify it
// under the terms found in the LICENSE[.md|.rst] file distributed
// together with this file.
//
// Lambda is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//
// ==========================================================================
// lambda_indexer.hpp: Main File for the indexer application
// ==========================================================================

#ifndef SEQAN_LAMBDA_LAMBDA_INDEXER_H_
#define SEQAN_LAMBDA_LAMBDA_INDEXER_H_

// #include <seqan/basic.h>
// #include <seqan/sequence.h>
//
// #include <seqan/seq_io.h>
// #include <seqan/index.h>
// #include <seqan/translation.h>
// #include <seqan/reduced_aminoacid.h>

#include "mkindex_misc.hpp"
#include "mkindex_saca.hpp"
// #include "shared_misc.hpp"
// #include "shared_options.hpp"
// #include "search_output.hpp" //TODO only needed because options are in one file, remove later

using namespace seqan;

inline void
printProgressBar(uint64_t & lastPercent, uint64_t curPerc)
{
    //round down to even
    curPerc = curPerc & ~1;
//     #pragma omp critical(stdout)
    if ((curPerc > lastPercent) && (curPerc <= 100))
    {
        for (uint64_t i = lastPercent + 2; i <= curPerc; i+=2)
        {
            if (i == 100)
                std::cout << "|" << std::flush;
            else if (i % 10 == 0)
                std::cout << ":" << std::flush;
            else
                std::cout << "." << std::flush;
        }
        lastPercent = curPerc;
    }
}

template <typename T>
inline void
myPrintImpl(//SharedOptions const & /**/,
            T const & first)
{
    std::cout << first;
}

inline void
myPrintImpl(//SharedOptions const & options,
            std::stringstream const & first)
{
    std::string str = first.str();
//     std::cerr << "terminal cols: " << options.terminalCols
//               << " str.size() " << str.size() << "\n";
    // if (options.isTerm && (str.size() >= (options.terminalCols -12)))
    //     std::cout << str.substr(str.size()-options.terminalCols+12,
    //                             options.terminalCols);
    // else
        std::cout << str;
}

template <typename T, typename ... Args>
inline void
myPrintImpl(//SharedOptions const & options,
            T const & first,
            Args const & ... args)
{
    myPrintImpl(/*options, */first);
    myPrintImpl(/*options, */args...);
}

template <typename ... Args>
inline void
myPrintImplThread(//SharedOptions const & options,
//                   T const & first,
                  Args const & ... args)
{
    SEQAN_OMP_PRAGMA(critical(stdout))
    {
//                 std::cout << "\033[" << omp_get_thread_num() << "B";
//                 std::cout << "\033E";
        // if (options.isTerm)
        // {
            for (unsigned char i=0; i< omp_get_thread_num(); ++i)
                std::cout << std::endl;
            std::cout << "\033[K";
        // }
        std::cout << "Thread " << std::setw(3) << omp_get_thread_num() << "| ";

        myPrintImpl(/*options, */args...);
        std::cout << "\n" << std::flush;
        // if (options.isTerm)
            std::cout << "\033[" << omp_get_thread_num()+1 << "A";
    }
}

template <typename... Args>
inline void
myPrint(/*SharedOptions const & options, const int verbose, */Args const &... args)
{
    // if (options.verbosity >= verbose)
    // {
        #if defined(_OPENMP)
        if (omp_in_parallel())
            myPrintImplThread(/*options,*/ args...);
        else
        #endif
            myPrintImpl(/*options,*/ args...);

        std::cout << std::flush;
    // }
}

// template <typename TRedAlph, BlastProgram p>
// void
// checkIndexSize(TCDStringSet<String<TRedAlph>> const & seqs,
//                LambdaIndexerOptions const & options,
//                BlastProgramSelector<p> const &)
// {
//     myPrint(options, 1, "Checking parameters of to-be-built index...");
//
//     // check number of sequences
//     using SAV = typename SAValue<TCDStringSet<String<TRedAlph>>>::Type;
//     uint64_t curNumSeq = length(seqs);
//     uint64_t maxNumSeq = std::numeric_limits<typename Value<SAV, 1>::Type>::max();
//
//     if (curNumSeq >= maxNumSeq)
//     {
//         throw std::invalid_argument(std::string("ERROR: Too many sequences to be indexed:\n  ") +
//                                     std::to_string(length(seqs)) +
//                                     std::string(" in file, but only ") +
//                                     std::to_string(maxNumSeq) +
//                                     std::string(" supported by index.\n"));
//     }
//
//     // check length of sequences
//     uint64_t maxLenSeq = std::numeric_limits<typename Value<SAV, 2>::Type>::max();
//     uint64_t maxLen = 0ul;
//     for (auto const & s : seqs)
//         if (length(s) > maxLen)
//             maxLen = length(s);
//
//     if (maxLen >= maxLenSeq)
//     {
//         std::string err;
//         err += "Sequences too long to be indexed:\n  ";
//         err += "length";
//         err += std::to_string(maxLen);
//         err += " present in file, but only ";
//         err += std::to_string(maxLenSeq);
//         err += " supported by index.\n";
//         #ifndef LAMBDA_LONG_PROTEIN_SUBJ_SEQS
//         if (p != BlastProgram::BLASTN)
//             err += "You can recompile Lambda and add -DLAMBDA_LONG_PROTEIN_SUBJ_SEQS=1 to activate\n"
//                    "support for longer protein sequences.\n";
//         #endif
//
//         throw std::invalid_argument(err);
//     }
//
//     // check available RAM
//     auto ram = getTotalSystemMemory();
//     auto lS = lengthSum(seqs);
//     unsigned long long factor = 0;
//     if (options.algo == "radixsort")
//         factor = sizeof(SizeTypeNum_<TRedAlph>) + sizeof(SizeTypePos_<TRedAlph>) + 4; // 4 is good heuristic
//     else if (options.algo == "skew7ext")
//         factor = 6; // TODO do some tests!
//     auto estimatedSize = lS * factor;
//
//     myPrint(options, 1, "done.\n");
//     if (estimatedSize >= ram)
//     {
//         std::cerr << "WARNING: Lambda estimates that it will need " << estimatedSize / 1024 / 1024 << "MB\n"
//                   << "         of memory to index this file, but you have only " << ram / 1024 / 1024 << "MB\n"
//                   << "         available on your system.\n"
//                   << "         This means you will likely encounter a crash with \"bad_alloc\".\n"
//                   << "         Split you sequence file into many smaller ones or use a computer\n"
//                   << "         with more memory!\n";
//     } else
//     {
//         myPrint(options, 2, "Detected RAM: ", ram / 1024 / 1024, "MB, Estimated RAM usage: ",
//                 estimatedSize / 1024 / 1024, "MB\n\n");
//     }
// }

// --------------------------------------------------------------------------
// Function createSuffixArray()
// --------------------------------------------------------------------------

// If there ids no overload with progress function, then strip it
template <typename TSA,
          typename TString,
          typename TSSetSpec,
          typename TAlgo,
          typename TLambda>
inline void
createSuffixArray(TSA & SA,
                  StringSet<TString, TSSetSpec> const & s,
                  TAlgo const &,
                  TLambda &&)
{
    return createSuffixArray(SA, s, TAlgo());
}

// ----------------------------------------------------------------------------
// Function indexCreate
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec, typename TConfig>
void
indexCreateProgress(Index<TText, FMIndex<TSpec, TConfig> > & index,
                    FibreSALF const &)
                    // LambdaIndexerOptions const & options)
{
    typedef Index<TText, FMIndex<TSpec, TConfig> >               TIndex;
    typedef typename Fibre<TIndex, FibreTempSA>::Type            TTempSA;
    typedef typename Size<TIndex>::Type                          TSize;
    typedef typename DefaultIndexCreator<TIndex, FibreSA>::Type  TAlgo;

    TText const & text = indexText(index);

    if (empty(text))
        return;

    TTempSA tempSA;
    uint64_t lastPercent = 0;

    double s = sysTime();
    myPrint(/*options, 1, */"Generating Index 0%  10%  20%  30%  40%  50%  60%  70%  80%  90%  100%\n"
                        " Progress:       |");
    // Create the full SA.
    resize(tempSA, lengthSum(text), Exact());
    // if (options.verbosity >= 1)
    // {
        createSuffixArray(tempSA,
                          text,
                          TAlgo(),
                          [&lastPercent] (uint64_t curPerc)
                          {
                              // needs locking, because called from multiple threads
                              SEQAN_OMP_PRAGMA(critical(progressBar))
                              printProgressBar(lastPercent, curPerc * 0.85); // 85% of progress
                          });
    // } else
    // {
    //     createSuffixArray(tempSA,
    //                       text,
    //                       TAlgo());
    // }
    double sacaTime = sysTime() - s;

    // if (options.verbosity >= 1)
        printProgressBar(lastPercent, 85);

    // Create the LF table.
    s = sysTime();
    // if (options.verbosity >= 1)
    // {
        createLFProgress(indexLF(index),
                         text,
                         tempSA,
                         [&lastPercent] (uint64_t curPerc)
                         {
                             // doesn't need locking, only writes from one thread
                             printProgressBar(lastPercent, curPerc * 0.1); // 10% of progress
                         });
    // } else
    // {
    //     createLFProgress(indexLF(index),
    //                      text,
    //                      tempSA,
    //                      [] (uint64_t) {});
    // }
    // Set the FMIndex LF as the CompressedSA LF.
    setFibre(indexSA(index), indexLF(index), FibreLF());
    double bwtTime = sysTime() - s;

    // if (options.verbosity >= 1)
        printProgressBar(lastPercent, 95);

    // Create the sampled SA.
    s = sysTime();
    TSize numSentinel = countSequences(text);
    createCompressedSa(indexSA(index), tempSA, numSentinel);
    double sampleTime = sysTime() - s;

    // if (options.verbosity >= 1)
        printProgressBar(lastPercent, 100);

    myPrint(/*options, 1,*/ "\n");
    myPrint(/*options, 2,*/ "SA  construction runtime: ", sacaTime, "s\n");
    myPrint(/*options, 2,*/ "BWT construction runtime: ", bwtTime, "s\n");
    myPrint(/*options, 2,*/ "SA  sampling runtime:     ", sampleTime, "s\n");
    myPrint(/*options, 1,*/ "\n");
}

template <typename TText, typename TSpec, typename TConfig>
void
indexCreateProgress(Index<TText, BidirectionalIndex<FMIndex<TSpec, TConfig> > > & index,
                    FibreSALF const &)
                    // LambdaIndexerOptions const & options)
{
    myPrint(/*options, 1,*/ "Bi-Directional Index [forward]\n");
    indexCreateProgress(index.fwd, FibreSALF()/*, options*/);

    myPrint(/*options, 1,*/ "Bi-Directional Index [backward]\n");
    indexCreateProgress(index.rev, FibreSALF()/*, options*/);
}

template <typename T>
inline void
_clearSparseSuffixArray(T &, std::false_type const &)
{}

template <typename T>
inline void
_clearSparseSuffixArray(T & dbIndex, std::true_type const &)
{
    // reverse index does not require sampled suffix array, but its size :|
    clear(getFibre(getFibre(getFibre(dbIndex, FibreSA()), FibreSparseString()), FibreValues()));
    clear(getFibre(getFibre(getFibre(dbIndex, FibreSA()), FibreSparseString()), FibreIndicators()));
}

#endif // header guard
