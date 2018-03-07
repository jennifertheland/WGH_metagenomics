#!/bin/bash

#
# Initials
#

processors=1
mailing=false
remove=false

#
# Argument parsing
#

while getopts "hp:" OPTION
do
    case ${OPTION} in

        p)
            processors=${OPTARG}
            ;;
        h)
            echo ''
	    echo 'This script performs fastq file processing and assembly of paired read metagenomics data.'
	    echo ""
	    echo 'Useage: bash athena_assembly.sh <r1_trimmed.fq> <r2_trimmed.fq> <NNN.clstr> <output_dir>'
	    echo ''
	    echo "Positional arguments (required)"
	    echo "  <r1_trimmed.fq>         Read one in .fastq format. The reads shold be trimmed with the barcode in the header. "
	    echo "  <r2_trimmed.fq>         Read two in .fastq format. The reads should be trimmed with the barcode in the header."
	    echo "  <NNN.clstr>     Concatenated .clstr file with all barcode clusters."
	    echo "  <output_dir>    Output directory holding all assembly files."
	    echo ""
	    echo "Optional arguments"
	    echo "  -h  help (this output)"
	    echo "  -p  processors for threading, not implemented yet"
	    echo ''
	    exit 0
	    ;;
    esac
done

ARG1=${@:$OPTIND:1} # r1 file
ARG2=${@:$OPTIND+1:1} # r2 file

file1=$ARG1
name_ext1=$(basename "$file1")
echo $name_ext1
name1="${name_ext1%.*}"
#file_name1="$path/${name_ext1%.*}"

file2=$ARG2
name_ext2=$(basename "$file2")
name2="${name_ext2%.*}" # name without extension

#
# Error handling
#

if [ -z "$ARG1" ] || [ -z "$ARG2" ]
then
    echo ""
    echo "ARGUMENT ERROR"
    echo "Did not find all two positional arguments, see -h for more information."
    echo ""
    exit 0
fi

dir_of_scripts=$(dirname "$0")
dir_of_files=$PWD


# head -n 40 in.fastq | awk '{if(NR%4==0) printf("%s",$0);}' | od -A n -t u1 | awk 'BEGIN{min=100;max=0;}{for(i=1;i<=NF;i++) {if($i>max) max=$i; if($i<min) min=$i;}}END{if(max<=74 && min<59) print "Phred+33"; else if(max>73 && min>=64) print "Phred+64"; else if(min>=59 && min<64 && max>73) print "Solexa+64"; else print "Unknown score encoding\!";}'


#for k in seq 20 8 124; do mkdir k$k abyss-pe -C $dir_of_files/k$k name=seed_opt k=$k in=$ARG1 $ARG2 done