## ARKS+LINKS pipeline example
## ARKS, LINKS and arks-make file must be in your path

if [ $# -ne 2 ]; then
    echo "Usage: $(basename $0) <draft> <reads>" >&2
    echo ""
    echo "draft		Assembled Sequences to further scaffold (Multi-Fasta format, with extension .fa or .fasta)"
    echo "reads		Interleaved chromium reads with barcode in the read header (ie. @ReadName BX:Z:<barcode>)"
    echo "		NOTE: file must have .fastq.gz or .fq.gz extension" 
    echo "This script expects the above files to be in the current working directory - Use soft links if needed."
    exit 1
fi

draft=$1; shift
reads=$1; shift

draft_base=$(basename $draft)
reads_base=$(basename $reads)

export draft_base
export reads_base

draft_prefix=$(perl -e 'if($ENV{'draft_base'} =~ /(\S+)\.fa/){print "$1";}')
reads_base=$(perl -e 'if($ENV{'reads_base'} =~ /(\S+)\.f(ast)?q\.gz/){print "$1";}')

echo $draft_prefix
echo $reads_base


##Running the arks makefile with default parameters (except specify 8 threads, -j 0.5 for ARKS):
arks-make  arks draft=${draft_prefix} reads=${reads_base} j=0.5 threads=8

