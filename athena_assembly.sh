#!/usr/bin/env bash

while getopts "hp:" OPTION
do
    case ${OPTION} in
        h)
            echo ''
	    echo 'This script processes paried read metagenomics data and performs de novo whole genome assembly..'
	    echo ""
	    echo 'Usage: bash athena_assembly.sh <r1_trimmed.fq> <r2_trimmed.fq> <NNN.clstr> <output_dir>'
	    echo ''
	    echo "Positional arguments (required)"
	    echo "  <r1_trimmed.fq>         Read one in .fastq format. The reads shold be trimmed with the barcode in the header. "
	    echo "  <r2_trimmed.fq>         Read two in .fastq format. The reads should be trimmed with the barcode in the header."
	    echo "  <NNN.clstr>     Concatenated .clstr file with all barcode clusters."
	    echo "  <output_dir>    Output directory holding all assembly files."
	    echo ""
	    echo "Optional arguments"
	    echo "  -h  help (this output)"
	    echo ''
	    exit 0
	    ;;
    esac
done

ARG1=${@:$OPTIND:1} # r1 file
ARG2=${@:$OPTIND+1:1} # r2 file
ARG3=${@:$OPTIND+2:1} # .clstr file
ARG4=${@:$OPTIND+3:1} # output dir

file1=$ARG1
name_ext1=$(basename "$file1")
name1="${name_ext1%.*}"

file2=$ARG2
name_ext2=$(basename "$file2")
name2="${name_ext2%.*}" # name without extension

if [ -z "$ARG1" ] || [ -z "$ARG2" ] || [ -z "$ARG3" ] || [ -z "$ARG4" ]
then
    echo ""
    echo "ARGUMENT ERROR"
    echo "Did not find all four positional arguments, see -h for more information."
    echo ""
    exit 0
fi

dir_of_scripts=$(dirname "$0")
dir_of_files=$PWD

source activate athena_assembly

#### STEP 1 ####

# TAG THE .FASTQ FILES WITH "BC:Z:BC_SEQ" USING tag_fastq.py

printf "### STEP 1 - Tagging initiated \n"

for file in $ARG1 $ARG2;
    do
    name_ext=$(basename "$file")
    printf 'Tagging %s \n' "$name_ext"
    name="${name_ext%.*}"
    python3 $dir_of_scripts/python_scripts/tag_fastq.py $file $ARG3 $name".tagged.fastq";
    done;

printf '### STEP 1 complete. Read files tagged with BC:Z:BC_SEQ \n'

#### STEP 2 ####

# SORT THE FILE ACCORDING TO BC_SEQ

printf "### STEP2 - Sorting initiated \n"

for file in $dir_of_files/*.tagged.fastq
    do
    name_ext=$(basename "$file")
    printf "Sorting %s \n" "$name_ext"
    name="${name_ext%.*}"
    bash $dir_of_scripts/sort_file.sh $file $name".sorted.fastq";
    done;

printf '### STEP 2 complete. Tagged read files sorted according to barcode sequence \n'

#### STEP 3 ####

# MERGE THE TAGGED AND SORTED .FASTQ FILES (CREATING AN INTERLEAVED FILE) AND CONVERT IT TO .FASTA FORMAT (bbmap reformat)

reformat.sh in1=$dir_of_files/$name1.tagged.sorted.fastq in2=$dir_of_files/$name2.tagged.sorted.fastq out=$dir_of_files/interleaved_R1_R2.fastq
reformat.sh in=$dir_of_files/interleaved_R1_R2.fastq out=$dir_of_files/interleaved_R1_R2.fasta

rm $dir_of_files/*tagged*

printf '### STEP 3 complete. read1.fastq and read2.fastq merged to interleaved file and converted to .fasta format. \n'

#### STEP 4 ####

# RUN IDBA TO ASSEMBLE SEED CONTIGS
printf "### STEP 4 - Initiating assembly of seed contigs \n"

mkdir idba_seed_contigs

idba_ud -r $dir_of_files/interleaved_R1_R2.fasta -o idba_seed_contigs

cp $dir_of_files/idba_seed_contigs/contig.fa $dir_of_files
rm -rf $dir_of_files/idba_seed_contigs
printf '### STEP 4 complete. Seed contigs generated.\n'

#### STEP 5 ####

# RUN BWA MEM TO MAP THE READ CLOUDS TO THE ASSEMBLED CONTIGS
printf "### STEP 5 - Initiating BWA alignment\n"
bwa index $dir_of_files/contig.fa
samtools faidx $dir_of_files/contig.fa

printf 'STEP 5.1 complete. BWA indexes generated.\n'

bwa mem -C -p $dir_of_files/contig.fa $dir_of_files/interleaved_R1_R2.fastq | \
    samtools sort -o mapped_reads.idba_contigs.bam -

printf 'STEP 5.2 complete. Reads mapped to seed contigs.\n'
printf '### STEP 5 complete.\n'

#### STEP 6 ####

# must make index of the .bam file in order for it to work.
printf "STEP 6 - Generating .bai file \n"
samtools index mapped_reads.idba_contigs.bam mapped_reads.idba_contigs.bam.bai
printf "STEP 6 complete. .bai file generated \n"
#### STEP 7 ####
# GENERATE CONFIG.JSON FILE

printf "STEP 7 - Generating config.json file \n"

echo -e "{" >> config.json
echo -e	"\t \"ctgfasta_path\" : \t \"$dir_of_files/contig.fa\"," >> config.json
echo -e	"\t \"reads_ctg_bam_path\":  \"$dir_of_files/mapped_reads.idba_contigs.bam\","  >> config.json
echo -e	"\t \"input_fqs\":            \"$dir_of_files/interleaved_R1_R2.fastq\"," >> config.json
echo -e	"\t \"cluster_settings\":  {" >> config.json
echo -e "\t\t \"cluster_type\": \"multiprocessing\"," >> config.json
echo -e	"\t\t \"processes\":" 4 >> config.json
echo -e	"\t}" >> config.json
echo -e  "}" >> config.json

printf "STEP 7 complete. \n"

#### STEP 8 ####
# RUN ATHENA

printf "STEP 8 - Running Athena-meta"

athena-meta $dir_of_files/config.json

printf "STEP 8 complete. Genome assembled."


#### STEP 9 ####

#bwa index $dir_of_files/results/olc/athena.asm.fa
#bwa mem -C -p $dir_of_files/results/olc/athena.asm.fa $dir_of_files/interleaved_R1_R2.fastq | \
#    samtools sort -n -o into_arcs.bam -

# echo "into_arcs.bam" >> bamfiles.txt

# arcs -f $dir_of_files/results/olc/athena.asm.fa -a bamfiles.txt -s 98 -c 5 -l 0 -z 500 -m 50-1000 -d 0 -e 1000 -r 0.05


# graph=athena.asm.fa.scaff_s98_c5_l0_d0_e1000_r0.05_original.gv

# python $dir_of_scripts/python_scripts/makeTSVfile.py $graph athena.asm.fa.c5_e1000_r0.05.tigpair_checkpoint.tsv $dir_of_files/results/olc/athena.asm.fa

