#!/usr/bin/env bash

processors=1

while getopts "hp:" OPTION
do
    case ${OPTION} in
        p)
           processors=${OPTARG}
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
	    echo '  -p <processors>     DEFAULT: 1'
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

path=$(dirname "$0") &&

printf '\n0. Argparsing & options' &&
printf '\nRead 1:\t'$ARG1'\nRead 2:\t'$ARG2'\nBCs:\t'$ARG3'\nOutput:\t'$ARG4'\n' &&


# START OF PROCESS #
printf 'BLR_metagenomics starting\n' &&
mkdir -p $output &&
#cd $output &&

### 1. READ PROCESSING ###
printf "### STEP 1 - Read processing \n" &&

# tag fastq files with clustered bc information
python3 $path/python_scripts/tag_fastq.py \
    $ARG1 \
    $ARG2 \
    $ARG3 \
    $file_name1'.tag.fastq' \
    $file_name2'.tag.fastq' &&

# sort and convert to itlvd fq
python3 $path/python_scripts/sort_tagged_fastq.py \
    -io \
    $file_name1'.tag.fastq' \
    $file_name2'.tag.fastq' \
    $file_name1'.tag.sort.itlvd.fastq' \
    'redundant_argument' &&

# fq => fa conversion
cat $file_name1'.tag.sort.itlvd.fastq' | paste - - - - | cut -f 1,2 | tr "\t" "\n" > $file_name1'.tag.sort.itlvd.fasta' &&

printf '### STEP 1 - Read processing complete\n' &&
### 2. READ CLOUD PRE-PROCESSING (make contig, map reads to contigs) ###
printf '### STEP 2 - Read cloud pre-processing\n' &&

idba_ud \
    -r $file_name1'.tag.sort.itlvd.fasta' \
    $file_name1'.preprocessing_assembly' &&

# IDBA automatically writes to a folder named out
assembly_output=$output'/mock_assembly' &&
mv 'out' $assembly_output &&

bwa index \
    $assembly_output'/contig.fa' &&

bwa mem \
    -C \
    -p $assembly_output'/contig.fa' \
    $file_name1'.tag.sort.itlvd.fastq' | \
    samtools sort -o $file_name1'.read_cloud_preprocessings.bam' - &&

printf '### STEP 2 - Read cloud pre-processing complete\n' &&
### 3. CONTIG GENERATION (Athena) ###
printf '### STEP 3 - Assembly of read clouds\n' &&

samtools index $file_name1'.read_cloud_preprocessings.bam' &&

printf  \
'{
    "ctgfasta_path" : "'$assembly_output'/contig.fa",
    "reads_ctg_bam_path" : "'$file_name1'.read_cloud_preprocessings.bam",
    "input_fqs" : "'$file_name1'.tag.sort.itlvd.fasta",

    "cluster_settings": {
        "cluster_type": "multiprocessing",
        "processes": '$processors'
    }
}' > $output'/config.json' &&

athena-meta $output'/config.json'



# scripts

#printf '### STEP 3 - Assembly of read clouds complete\n' &&
### 4. SCAFFOLDING ###
#printf '### STEP 4 - Contig scaffolding\n' &&
# scripts

#printf '### STEP 4 - Contig scaffolding complete\n' &&

# POST ANALYSIS PRINTING ETC.
#printf 'BLR_metagenomics FINISHED\n' &&

#cd -











#printf "### STEP 4 - Initiating assembly of seed contigs \n"
#mkdir idba_seed_contigs
#idba_ud -r $dir_of_files/interleaved_R1_R2.fasta -o idba_seed_contigs
#
#cp $dir_of_files/idba_seed_contigs/contig.fa $dir_of_files
#rm -rf $dir_of_files/idba_seed_contigs
#
#printf '### STEP 4 complete. Seed contigs generated.\n'
#
## RUN BWA MEM TO MAP THE READ CLOUDS TO THE ASSEMBLED CONTIGS
#printf "### STEP 5 - Initiating BWA alignment\n"
#bwa index $dir_of_files/contig.fa
#samtools faidx $dir_of_files/contig.fa
#
#printf 'STEP 5.1 complete. BWA indexes generated.\n'
#
#bwa mem -C -p $dir_of_files/contig.fa $dir_of_files/interleaved_R1_R2.fastq | \
#    samtools sort -o mapped_reads.idba_contigs.bam -
#
#printf 'STEP 5.2 complete. Reads mapped to seed contigs.\n'
#printf '### STEP 5 complete.\n'
#
#
#
## must make index of the .bam file in order for it to work.
#printf "STEP 6 - Generating .bai file \n"
#samtools index mapped_reads.idba_contigs.bam mapped_reads.idba_contigs.bam.bai
#printf "STEP 6 complete. .bai file generated \n"
##### STEP 7 ####
## GENERATE CONFIG.JSON FILE
#
#printf "STEP 7 - Generating config.json file \n"
#
#echo -e "{" >> config.json
#echo -e	"\t \"ctgfasta_path\" : \t \"$dir_of_files/contig.fa\"," >> config.json
#echo -e	"\t \"reads_ctg_bam_path\":  \"$dir_of_files/mapped_reads.idba_contigs.bam\","  >> config.json
#echo -e	"\t \"input_fqs\":            \"$dir_of_files/interleaved_R1_R2.fastq\"," >> config.json
#echo -e	"\t \"cluster_settings\":  {" >> config.json
#echo -e "\t\t \"cluster_type\": \"multiprocessing\"," >> config.json
#echo -e	"\t\t \"processes\":" 4 >> config.json
#echo -e	"\t}" >> config.json
#echo -e  "}" >> config.json
#
#printf "STEP 7 complete. \n"
#
##### STEP 8 ####
#
## RUN ATHENA
#
#printf "STEP 8 - Running Athena-meta \n"
#
#athena-meta $dir_of_files/config.json
#
#printf "STEP 8 complete. Genome assembled.\n"
#
#
##### STEP 9 ####
#
## RUN ARCS AND LINKS
#
#printf "STEP 9.1 - Aligning linked-reads to the assembled genome \n"
#
#bwa index $dir_of_files/results/olc/athena.asm.fa
#bwa mem -C -p $dir_of_files/results/olc/athena.asm.fa $dir_of_files/interleaved_R1_R2.fastq | \
#    samtools sort -n -o into_arcs.bam -
#
#echo "into_arcs.bam" >> bamfiles.txt
#
#printf "STEP 9.1 complete. \n"
#
#printf "STEP 9.2 - Running arcs \n"
#
#arcs -f $dir_of_files/results/olc/athena.asm.fa -a bamfiles.txt -s 98 -c 5 -l 0 -z 500 -m 50-1000 -d 0 -e 5000 -r 0.05
#
#graph=athena.asm.fa.scaff_s98_c5_l0_d0_e1000_r0.05_original.gv
#
#python $dir_of_scripts/python_scripts/makeTSVfile.py $graph athena.asm.fa.c5_e1000_r0.05.tigpair_checkpoint.tsv $dir_of_files/results/olc/athena.asm.fa
#
#printf "STEP 9.2 complete.\n"
#
#printf "STEP 9.3 - Running LINKS \n"
#
#touch empty.fof
#
#LINKS -f $dir_of_files/results/olc/athena.asm.fa -s empty.fof -k 20 -b athena.asm.fa.c5_e30000_r0.05 -l 5 -t 2 -a 0.3
#
#printf "STEP 9.3 complete. \n"
#printf "STEP 9 complete. Scaffolding done.\n"