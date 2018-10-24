# Barcode-Linked Reads Analysis

This GitHub repository describes and distributes the metagenomics pipeline used in "[Efficient whole genome haplotyping and 
high-throughput single molecule phasing with barcode-linked reads](https://www.biorxiv.org/content/early/2018/06/26/356121)".

## Dependencies

Here follows a list with links to all bioinformatics software needed to use of this part of the pipeline.

  - [IDBA-UD](https://github.com/loneknightpy/idba)
  - [Athena-meta](https://github.com/abishara/athena_meta)
  - [BWA](https://sourceforge.net/projects/bio-bwa/files/)
  - [LINKS](https://github.com/bcgsc/LINKS)
  - [arcs](https://github.com/bcgsc/arcs)
  - [samtools](https://github.com/samtools/samtools)
 
To utilize all aspects of the pipeline some GNU software are also needed.

  - [pigz](https://zlib.net/pigz/)
  - [mail](https://mailutils.org/manual/mailutils.html)

## Setup

First, download this GitHub repository by writing the cloning command in your terminal.

```
git clone https://github.com/jennifertheland/BLR_metagenomics.git
```

## Useage

The main script to run the pipeline is athena_assembly.sh, for available options, see -h (--help) and for more 
detail of what is run consult the step-by-step folder.

### Pre-processing

First, make sure you have trimmed your reads and clustered your barcodes (using the [BLR pipeline](https://github.com/FrickTobias/BLR)).

```
bash BLR_automation.sh -e 2 -m john.doe@myworkplace.com -p <processors> -r <read1.fq> <read2.fq> <read_processing_folder>
```

This will yield you two trimmed read files and one clustered barcode file which will be used as input into 
the metagenomics pipeline.

```
ls -lh <read_processing_folder>
1_trim.log
1_cluster.log
read1.trimmed.fq.gz
read2.trimmed.fq.gz
BC.clstr
```

### Assembly and scaffolding

To run the pipeline, write the following to yield an output folder with the results in.

```
bash athena_assembly.sh -m john.doe@myworkplace.com -p <processors> -m  <read1.trimmed.fq.gz> <read2.trimmed.fq.gz> <BC.NNN.clstr> <output>
```