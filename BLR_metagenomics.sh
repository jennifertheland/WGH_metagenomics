#!/usr/bin/env bash
set -euo pipefail

processors=1
single_lib=true
only_tag=false

while getopts "hp:mo" OPTION
do
    case ${OPTION} in
        p)
           processors=${OPTARG}
           ;;
        o)

           only_tag=true
           ;;
        m)
           single_lib=false
           ;;
        h)
        echo ''
	    echo 'This script processes paried read metagenomics data and performs de novo whole genome assembly..'
	    echo ""
	    echo 'Usage: bash BLR_metagenomics.sh <r1_trimmed.fq> <r2_trimmed.fq> <NNN.clstr> <output_dir>'
	    echo ''
	    echo "Positional arguments (required)"
	    echo "  <r1_trimmed.fq>         R1 in .fq format. Trimmed with the barcode in the header. "
	    echo "  <r2_trimmed.fq>         Read two in .fastq format. The reads should be trimmed with the barcode in the header."
	    echo "  <NNN.clstr>     Concatenated .clstr file with all barcode clusters."
	    echo "  <output_dir>    Output directory holding all assembly files."
	    echo ""
	    echo "Optional arguments"
	    echo "  -h  help (this output)"
	    echo '  -p  <processors>        DEFAULT: 1'
	    echo '  -o  omit tagging        DEFAULT: false'
	    echo '  -m  multiple libraries  DEFAULT: false'
	    exit 0
	       ;;
    esac
done

ARG1=${@:$OPTIND:1} # r1 file
ARG2=${@:$OPTIND+1:1} # r2 file
ARG3=${@:$OPTIND+2:1} # .clstr file
ARG4=${@:$OPTIND+3:1} # output dir

output=$ARG4
mkdir -p $output

file1=$ARG1
name_ext1=$(basename "$file1")
name1="${name_ext1%.*}"
file_name1="$output/${name_ext1%.*}"

file2=$ARG2
name_ext2=$(basename "$file2")
name2="${name_ext2%.*}" # name without extension
file_name2="$output/${name_ext2%.*}"


file3=$ARG3
name_ext3=$(basename "$file3")
name3="${name_ext3%.*}" # name without extension
file_name3="$output/${name_ext3%.*}"

if [ -z "$ARG1" ] || [ -z "$ARG2" ] || [ -z "$ARG3" ] || [ -z "$ARG4" ]
then
    echo ""
    echo "ARGUMENT ERROR"
    echo "Did not find all four positional arguments, see -h for more information."
    echo "(got r1:"$ARG1", r2:"$ARG2", BCs:"$ARG3" and output:"$ARG4" instead)"
    echo ""
    exit 0
fi

path=$(dirname "$0")

printf '\n0. Argparsing & options'
printf '\nRead 1:\t'$ARG1'\nRead 2:\t'$ARG2'\nBCs:\t'$ARG3'\nOutput:\t'$ARG4'\n'

#
# # # END OF PARSING & SETUP
#

# START OF PROCESS #
printf '\nBLR_metagenomics starting\n'
mkdir -p $output

### 1. READ PROCESSING ###

# If multiple libraries are run, don't tag (should be done prior, manually)
if $single_lib
then
    printf "\n### STEP 1 - Read processing \n"
    # tag fastq files with clustered bc information
    python3 $path/python_scripts/tag_fastq.py \
        $ARG1 \
        $ARG2 \
        $ARG3 \
        $file_name1'.tag.fastq' \
        $file_name2'.tag.fastq'
else
    echo ''
    echo '!!! UNSATBLE SOLUTION: Copying inputs so pipeline finds the pre-tagged read files.'
    echo ''
    cp $ARG1 $file_name1'.tag.fastq'
    cp $ARG2 $file_name2'.tag.fastq'
    printf "### STEP 1 - Read processing \n"
    printf "Multiple library option active - omitting tagging step (should be done prior with -o option) \n"
fi

# If running multiple libraries, should stop here after tagging to combine libraries.
if $only_tag
then
    printf "\n### ONLY TAG OPTION ACTIVE - EXITING ###\n\n"
    exit
fi

# sort and convert to itlvd fq
python3 $path/python_scripts/sort_tagged_fastq.py \
    -io \
    $file_name1'.tag.fastq' \
    $file_name2'.tag.fastq' \
    $file_name1'.tag.sort.itlvd.fastq' \
    'redundant_argument'

# fq => fa conversion
cat $file_name1'.tag.sort.itlvd.fastq' | paste - - - - | cut -f 1,2 | tr "\t" "\n" > $file_name1'.tag.sort.itlvd.fasta'

printf '### STEP 1 - Read processing complete\n'
### 2. READ CLOUD PRE-PROCESSING (make contig, map reads to contigs) ###
printf '### STEP 2 - Read cloud pre-processing\n'

idba idba_ud \
    -r $file_name1'.tag.sort.itlvd.fasta' \
    $file_name1'.preprocessing_assembly'

# IDBA automatically writes to a folder named out
assembly_output=$output'/mock_assembly'
mv 'out' $assembly_output

bwa index \
    $assembly_output'/contig.fa'

bwa mem \
    -C \
    -p $assembly_output'/contig.fa' \
    $file_name1'.tag.sort.itlvd.fastq' | \
    samtools sort -o $file_name1'.read_cloud_preprocessings.bam' -

printf '### STEP 2 - Read cloud pre-processing complete\n'
### 3. CONTIG GENERATION (Athena) ###
printf '### STEP 3 - Assembly of read clouds\n'

samtools index $file_name1'.read_cloud_preprocessings.bam'

printf  \
'{
    "ctgfasta_path" : "'$assembly_output'/contig.fa",
    "reads_ctg_bam_path" : "'$file_name1'.read_cloud_preprocessings.bam",
    "input_fqs" : "'$file_name1'.tag.sort.itlvd.fasta",

    "cluster_settings": {
        "cluster_type": "multiprocessing",
        "processes": '$processors'
    }
}' > $output'/config.json'

athena-meta $output'/config.json'

printf '### STEP 3 - Assembly of read clouds done\n'
# Run flye
printf '### STEP 4 - Flye, NOT IMPLEMENTED INTO PIPELINE YET\n'