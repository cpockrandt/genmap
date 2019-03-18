#!/bin/sh

#echo "Downloading an indexed version of GRCh38 (5.0 GB) with dna4 alphabet ..."
#wget http://ftp.imp.fu-berlin.de/pub/cpockrandt/grch38-dna4.tar.gz

#echo "Uncompressing grch38-dna4.tar.gz ..."
#tar -zxvf grch38-dna4.tar.gz

IFS=','

if [ "$#" -ne 4 ]; then
    echo "USAGE:"
    echo "./bench.sh BIN_TO_GENMAP INDEX_DIRECTORY OUTPUT_DIRECTORY THREADS"
    exit 1
fi

BIN=$1
INDEXDIR=$2
OUTPUTDIR=$3
THREADS=$4

STATS=""

runGenMap() {
    K=$1; E=$2
    echo "($K, $E)-frequency:"
    /usr/bin/time -f "%e" -o $OUTPUTDIR/timing $BIN map -I $INDEXDIR -O $OUTPUTDIR -E $E -K $K -fs -r -m -T $THREADS 
    tput cuu 3
    echo -n "($K, $E)-frequency: "
    cat $OUTPUTDIR/timing
    #echo "----------------------"
    #STATS+="($K, $E)-frequency: 0.00s\n"
}

for i in 5,0 6,0; do # 36,0 24,1 36,2 50,2 75,3; do
    set -- $i
    runGenMap $1 $2
done

for i in 101,0 101,1 101,2 101,3 101,4; do
    set -- $i
    runGenMap $1 $2
done
