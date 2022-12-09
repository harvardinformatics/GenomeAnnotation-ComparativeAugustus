# AUGUSTUS: multi-species mode
For years, AUGUSTUS has been a key component of genome annotation pipelines, and has also been used as a stand-alone genome annotation tool. This page describes how to implement multi-genome, "comparative" AUGUSTUS (CGP). CGP uses a whole genome alignment (WGA) in order to leverage extrinsic evidence and for jointly making gene predictions across the genomes in the alignment. As CGP uses a whole-genome alignment as its backbone, it will not annotate genomic intervals in the reference genome specified in the MAF file for which there is no alignment. In other words, if a species of interest has functional gene sequence that doesn't align to the reference, it will be missed. Therefor, the AUGUSTUS developers recommending integrating outputs from single-genome and multi-species instances of AUGUSTUS to maximize sensitivity. Single-genome implementations of AUGUSTUS are now optimally implemented in [BRAKER](https://github.com/Gaius-Augustus/BRAKER), and we provide guidelines for running BRAKER in an HPC environment [here](https://github.com/harvardinformatics/GenomeAnnotation-Braker). It should be noted that there are myriad command line options for running CGP, and some not included here may be relevant depending upon the nature of the study being considered. 

## Converting a Cactus hal to MAF
The instructions on this page for running CGP are based upon whole-genome alignments obtained with [Cactus](https://github.com/ComparativeGenomicsToolkit/cactus), [*Armstrong et al. 2020*](https://www.nature.com/articles/s41586-020-2871-y), and require converting a Cactus [hal](https://github.com/ComparativeGenomicsToolkit/hal) binary to a multiple alignment format [MAF](https://genome.ucsc.edu/FAQ/FAQformat.html#format5). This can be done following [recommendations](https://bioinf.uni-greifswald.de/augustus/binaries/tutorial-cgp/cactus.html#hal2maf) provided by the AUGUSTUS developers. In principle, WGA produced by other bioinformatics tools and in formats other than hal can also be used. A snakemake workflow for generating a WGA with GPU-enabled Cactus can be found [here](https://github.com/harvardinformatics/GenomeAnnotation-WholeGenomeAlignment).

## CGP with protein evidence 
### Hint creation
To generate splice hints from external protein sequence data, we use scripts that are part of the Augustus distribution to wrap protein alignment to the genome using [GenomeThreader](https://genomethreader.org/) and creation of hints from those alignments. 

We first create a genomethreader conda environment
```bash
module load python
conda create -n genomethreader -c bioconda genomethreader
```

Next we execute the Augustus startAlign.pl script:
```bash
#!/bin/sh
module load python
source activate genomethreader

# be sure to specify the full path to the perl script
genomefasta=$1
proteintargetfasta=$2

startAlign.pl --genome=${genomefasta} --prot=${proteintargetfasta} --CPU=12 --prg=gth
```


where $genomefasta and $proteintargetfasta are the genome sequences for the species of interest, and the protein fasta from (ideally) closely related species. For a large protein fasta, one can split the fasta file into smaller pieces, run aligngment on those pieces, and combine the alignment files afterwards before proceeding to the next step. The number of CPUs flag should be changed to reflect the available compute resources you have at your disposal.

We run the perl script from outside of the singularity container we use for running AUGUSTUS to avoid various path-related errors. This analysis will produce a directory with the name align_gth. Next, we convert the alignment to splice hints. To do this, download the BRAKER [align2hints.pl](https://github.com/Gaius-Augustus/BRAKER/blob/master/scripts/align2hints.pl) script. Then run it as follows:
```bash
align2hints.pl --in align_gth/gth.concat.aln --out=prot.hints --prg=gth
```
This script looks for the align_gth directory produced during the alignment step so run it from the same directory in which align_gth is located. The file *gth.concat.aln* is the file containing all of the GenomeThreader alignments, i.e. if you ran alignments in parallel on many subsets of the protein fasta, they should be concatenated to produce *gth.concat.aln*. 

## CGP with RNA-seq evidence
### Hint creation
To generate RNA-seq derived splice-site hints for Augustus is a multi-step process which present computational challenges. We have implemented changes to tools within the Augustus distribution in order for them to handle contemporary, increasingly large RNA-seq datasets in a time-efficient manner. 

After using samtools to merge the bamfiles of interest, one has to sort them bam-file by sequence name rather that coordinate order. In the directory where the merged bam is, first create a tmp directory:
```bash
mkdir tmp
```
Next, use samtools to namesort the bam file, using our script [namesortbam.sh]():
```bash
sbatch namesortbam.sh $mymergedbam
```

After name-sorting, the bam file needs to be filtered. We have modified the original [filterBam](https://github.com/nextgenusfs/augustus/tree/master/auxprogs/filterBam) code to [filterBam-zlib-ng-2:LINK NOT YET ADDED] in order to speed up filtering. We execute this updated version from within a singularity image,oneapi-hpckit_2021.2-devel-centos8.sif, available at: . To execute this step we use [filterbam_FAS_informatics_version.sh](https://github.com/harvardinformatics/GenomeAnnotation/blob/master/ComparativeAugustus/slurm_scripts/filterbam_FAS_informatics_version.sh) run as follows:
```bash
sbatch filterbam_FAS_informatics_version.sh $mymergeddbam
```

After filtering the bam file, we then generate "intron part" hints by running the Augustus tool bam2hints with [bam2hints](https://github.com/harvardinformatics/GenomeAnnotation/blob/master/ComparativeAugustus/slurm_scripts/bam2hints.sh):
```bash
sbatch bam2hints $myfilteredbam
``` 

To generate "exon part" hints, we then resort the filtered bam file back into coordinate order with [coordsortbam.sh](https://github.com/harvardinformatics/GenomeAnnotation/blob/master/ComparativeAugustus/slurm_scripts/coordsortbam.sh):
```bash
mkdir tmp
sbatch coordsortbam.sh $myfilteredbam
```

Then, we convert the resorted bam to wig format with the Augustus tool bam2wig using [bam2wig.sh](https://github.com/harvardinformatics/GenomeAnnotation/blob/master/ComparativeAugustus/slurm_scripts/bam2wig.sh):
```bash
sbatch bam2wig.sh $coordsorted_filtered_bam
```
Finally, we convert the wig file into a gff formatted hints file with wig2hints.pl via our slurm script [wig2hints.sh](https://github.com/harvardinformatics/GenomeAnnotation/blob/master/ComparativeAugustus/slurm_scripts/wig2hints.sh)
```bash
sbatch wig2hints.sh $mywigfile
```
