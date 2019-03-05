#!/bin/sh

errorout()
{
    echo $1 #> /dev/stderr
    [ "$MYTMP" = "" ] || rm -r "${MYTMP}"
    exit 1
}

[ $# -ne 6 ] && exit 1

SRCDIR=$1
BINDIR=$2
CASE=$3
INDEX_FLAGS=$4
FLAGS=$5
EXPECTED_FOLDER=$6

# check existence of commands
which mktemp diff > /dev/null
[ $? -eq 0 ] || errorout "Not all required programs found. Needs: mktemp diff"

MYTMP="$(mktemp -q -d -t "$(basename "$0").XXXXXX" 2>/dev/null || mktemp -q -d)"
[ $? -eq 0 ] || errorout "Could not create tmp"

mkdir -p "${MYTMP}/output"
[ $? -eq 0 ] || errorout "Could not create folder in tmp"

cd "$MYTMP"
[ $? -eq 0 ] || errorout "Could not cd to tmp"

if [ "$INDEX_FLAGS" = "-FD" ]; then
    ${BINDIR}/bin/genmap index -FD "${SRCDIR}/tests/test_cases/case_${CASE}" -I "${MYTMP}/index" -A skew
else
    ${BINDIR}/bin/genmap index -F "${SRCDIR}/tests/test_cases/case_${CASE}/genome.fa" -I "${MYTMP}/index" -A skew
fi

${BINDIR}/bin/genmap map -I "${MYTMP}/index" -O "${MYTMP}/output" ${FLAGS}
diff -r --strip-trailing-cr "${SRCDIR}/tests/test_cases/case_${CASE}/${EXPECTED_FOLDER}" "${MYTMP}/output"
[ $? -eq 0 ] || errorout "Files are not equal!"

${BINDIR}/bin/genmap map -I "${MYTMP}/index" -O "${MYTMP}/output" ${FLAGS} -xo 1
diff -r --strip-trailing-cr "${SRCDIR}/tests/test_cases/case_${CASE}/${EXPECTED_FOLDER}" "${MYTMP}/output"
[ $? -eq 0 ] || errorout "Files are not equal!"

${BINDIR}/bin/genmap map -I "${MYTMP}/index" -O "${MYTMP}/output" ${FLAGS} -xo 2
diff -r --strip-trailing-cr "${SRCDIR}/tests/test_cases/case_${CASE}/${EXPECTED_FOLDER}" "${MYTMP}/output"
[ $? -eq 0 ] || errorout "Files are not equal!"

# gunzip < "${SRCDIR}/tests/db_${SALPHIN}.fasta.gz" > db.fasta
# [ $? -eq 0 ] || errorout "Could not unzip database file"
#
# ${BINDIR}/bin/lambda_indexer -d db.fasta -di ${DI} -p ${PROG}
# [ $? -eq 0 ] || errorout "Could not run the indexer"
#
# openssl md5 * > md5sums
# [ $? -eq 0 ] || errorout "Could not run md5 or md5sums"
#
# gunzip < "${SRCDIR}/tests/db_${SALPH}_${DI}.md5sums.gz" > md5sums.orig
# [ $? -eq 0 ] || errorout "Could not unzip md5sums.orig"
#
# [ "$(cat md5sums)" = "$(cat md5sums.orig)" ] || errorout "$(diff -u md5sums md5sums.orig)"
#
# ## INDEXER tests end here
# if [ "$MODE" = "MKINDEX" ]; then
#     rm -r "${MYTMP}"
#     exit 0
# fi
#
# gunzip < "${SRCDIR}/tests/queries_${QALPHIN}.fasta.gz" > queries.fasta
# [ $? -eq 0 ] || errorout "Could not unzip queries.fasta"
#
# ${BINDIR}/bin/lambda -d db.fasta -di ${DI} -p ${PROG} -q queries.fasta -t 1 --version-to-outputfile off \
# -o output_${PROG}_${DI}.${EXTENSION}
# [ $? -eq 0 ] || errorout "Search failed."
#
# [ "$(openssl md5 output_${PROG}_${DI}.${EXTENSION})" = \
# "$(zgrep "(output_${PROG}_${DI}.${EXTENSION})" "${SRCDIR}/tests/search_test_outfile.md5sums.gz")" ] || errorout "MD5 mismatch of output file"

rm -r "${MYTMP}"
[ $? -eq 0 ] || errorout "Could not remove tmp"
