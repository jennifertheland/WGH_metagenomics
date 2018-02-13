#!/bin/bash

#
# Initials
#

processors=1

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
	    echo 'This script sorts a tagged .fastq file based on barcode sequence.'
	    echo ""
	    echo 'Usage: bash sort_file.sh <tagged.fastq> <tagged_sorted.fastq> '
	    echo ''
	    echo "Positional arguments (required)"
	    echo "  <tagged.fastq>         Read file in .fastq format. Should be tagged according to RG:Z:X BC:Z:X"
	    echo "  <tagged_sorted.fastq>  Sorted .fastq file, according to barcode sequence."
	    echo ""
	    echo "Optional arguments"
	    echo "  -h  help (this output)"
	    echo "  -p  processors for threading, not implemented yet"
	    echo ''
	    exit 0
	    ;;
    esac
done

#
# Positional redundancy for option usage
#

ARG1=${@:$OPTIND:1}
ARG2=${@:$OPTIND+1:1}

wgh_path=$(dirname "$0")

#
# Error handling
#

if [ -z "$ARG1" ] || [ -z "$ARG2" ]
then
    echo ""
    echo "ARGUMENT ERROR"
    echo "Did not find all positional arguments, see -h for more information."
    echo "(got input:"$ARG1" and output:"$ARG2" instead)"
    echo ""
    exit 0
fi

cat $ARG1 |

while read first_line; read second_line; read third_line; read fourth_line
    do
        echo "$first_line" "$second_line" "$third_line" "$fourth_line" >> temp.fastq;
    done

echo "### Read records collapsed into one line - DONE"

cat temp.fastq | sort -d -k2,2 -k1,1 > temp2.fastq

echo "### Reads sorted according to barcode sequence - DONE"

rm temp.fastq

echo "### Initiating reformatting to .fastq format."

cat temp2.fastq |

python $wgh_path/python_scripts/sort_tagged_file.py temp2.fastq $ARG2

echo "### Sorting - DONE"

rm temp2.fastq

