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

Also make sure your IDBA-UD installation can work on 150 bp read constructs by changing the kMaxShortSequence to
200 bases in the short_sequences.h file located in your idba folder under.

```
idba/src/sequence/short_sequence.h

```
Change line 102

```
static const uint32_t kMaxShortSequence = 128;

```
To this:

```
static const uint32_t kMaxShortSequence = 200;
```

After changing the line, remember to rebuild your IDBA installation.

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
bash BLR_metagenomics.sh -m john.doe@myworkplace.com -p <processors> <read1.trimmed.fq.gz> <read2.trimmed.fq.gz> <BC.NNN.clstr> <output>
```

## Advanced usage: Combine multiple libries

Start by running the pre-processing (see above) for all your libraries separately (the example shows two libraries but it 
can be done for any number of libraries).

```
# Library 1
bash BLR_automation.sh -e 2 -m john.doe@myworkplace.com -p <processors> -r <read1.fq> <read2.fq> <read_processing_folder>
# Library 2
bash BLR_automation.sh -e 2 -m john.doe@myworkplace.com -p <processors> -r <read1.fq> <read2.fq> <read_processing_folder>
```

Proceed by tagging all libraries with error corrected barcodes by using the -o (only tag) flag. 

```
# Library 1
bash BLR_metagenomics.sh -o <lib1.read1.trim.fq> <lib1.read2.trim.fq> <BC.NNN.clstr> <lib1_out>
# Library 2
bash BLR_metagenomics.sh -o <lib2.read1.trim.fq> <lib2.read2.trim.fq> <BC.NNN.clstr> <lib2_out>
```

Now combine the libraries using the multilib_combiner.py script found in the python_scripts folder.

```
python3 python_scripts/multilib_combiner.py <lib1_out/read1.tag.fq> <lib1_out/read2.tag.fq> <lib2_out/read1.tag.fq> <lib2_out/read2.tag.fq> -r1 <comb.read1.tag.fq> -r2 <comb.read2.tag.fq>
```

Lastly run the automation script with the -m (multiple libraries) flag.

```
bash BLR_metagenomics.sh -m <comb.read1.tag.fq> <comb.read2.tag.fq> <BC.NNN.clstr> <lib1_out>
```