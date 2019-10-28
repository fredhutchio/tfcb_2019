MCB517A Homework 3
================
Student Name
10/22/2019

Complete this homework by writing R code to complete the following
tasks. Keep in mind:

1.  Empty chunks have been included where code is required
2.  For Problem 2e, you should include a image (screen shot) instead of
    providing code
3.  This homework requires use of
    `BRCA.genome_wide_snp_6_broad_Level_3_scna.seg`. If you choose to
    include this file in the same directory as your homework file,
    please do not include it in your homework submission (e.g., do not
    commit to GitHub).
4.  You will be graded on your code and output results (knitted .html
    file). The assignment is worth 40 points total; partial credit can
    be awarded.

This assignment is due at 12pm on Oct 29, 2019. Your responses should be
submitted through the private repository created via GitHub Classroom
using the link emailed to you. See the note at the bottom of this
document for more information on formatting your final
submission.

## Problem 1: Overlaps between genomic regions and copy number alterations. (4 points total)

### Preparation

Load copy number segment results as shown in *2.1 BED format* of
*Lecture8\_GenomicData.Rmd*. You will use the same file as in the
lecture notes, `BRCA.genome_wide_snp_6_broad_Level_3_scna.seg`. Here is
code to get you started.

``` r
library(GenomicRanges)
library(data.table)
segs <- read.delim("BRCA.genome_wide_snp_6_broad_Level_3_scna.seg", as.is = TRUE)
mode(segs$Chromosome) <- "character"
segs[segs$Chromosome == 23, "Chromosome"] <- "X"
segs.gr <- as(segs, "GRanges")
```

### a. Find the segments in `segs.gr` that have *any* overlap with the region `chr8:128,746,347-128,755,810` (2 points)

Print out the first five unique TCGA IDs after sorting
alphabetically/lexicographically using the `sort` function with
`decreasing =
FALSE`.

### b. Find the mean of the `Segment_Mean` values for copy number segments that have *any* overlap with the region chr17:37,842,337-37,886,915. (2 points)

## Problem 2: Frequency of copy number alteration events within genomic regions. (12 points total)

This problem will continue to use the copy number data stored in
`segs.gr`.

### a. Create a genome-wide tile of 1Mb windows for chromosomes 1-22 and X. (2 points)

### b. Count the number of copy number alteration segments overlapping each 1Mb tile in chromosome X. (2 points)

For each 1Mb tile in chrX, print out the count of `any` type of overlap
of copy number segments in
`segs.gr`.

### c. Find the 1Mb window with the most *complete* overlap of segments (3 points).

Here, we want to find the overlap of segments from `segs.gr` that is
*completely* `within` a 1Mb bin created in (a). Return the 1Mb window
`GRanges` entry with the highest count of segments fully contained
`within`. Hints: see *3.3 Counting overlapping ranges* of
*Lecture8\_GenomicData.Rmd* if you want to use `data.table` for
this.

### d. Find the 1Mb window with the most frequent overlapping deletions. (3 points)

Find the 1Mb windows with `any` overlap with deletion copy number
segments. Assume a deletion segment is defined as a segment with
`Segment_Mean < -0.3`. Return the 1Mb window `Granges` entry with the
highest count of deletion
segments.

### e. Visually inspect the deletion overlap result from part (d) using IGV. (2 points)

Provide a screen shot of IGV at the 1Mb window with the most frequent
overlap with deletion segments. The image should include the segments
from `RCA.genome_wide_snp_6_broad_Level_3_scna.seg` loaded.

*Include image here*

## Problem 3: Finding Copy Number of Genes (8 points total)

This problem will continue to use the copy number data stored in
`segs.gr`.

### Preparation

You will use the Ensembl gene annotations from `biomaRt` to annotate
copy number segments. Here is some code to help you get started with
loading the `GRanges` for the `protein_coding` gene list.

``` r
library(biomaRt)
ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
annot <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol",
                          "entrezgene_id", "chromosome_name",
                          "start_position" ,"end_position",
                          "gene_biotype"),
             filters = "biotype", values = "protein_coding",
             mart = ensembl)
```

### a. Convert the data.frame `annot` into a sorted `GRanges` object with only chromosomes 1-22 and X. (2 points)

Hint: for sorting, use `sort` directly on the `GRanges`
object.

### b. Find the segments in `segs.gr` that have `any` overlap with the gene `ERBB2`. (2 points)

Print the first five `Granges` entries in `segs.gr`, sorted by range
coordinates, that overlap the gene coordinates of `ERBB2`. Hint: for
sorting, use `sort` directly on the `GRanges` object.

### c. Find the distribution of copy number for `PIK3CA`. (4 points)

Find the counts of deletion (`Segment_Mean < -0.3`), neutral
(`Segment_Mean >= -0.3 & Segment_Mean <= 0.3`), gain (`Segment_Mean> 0.3`) 
segments that have `any` overlap with `PIK3CA` gene coordinates.

## Problem 4: Reading and annotating genomic variants (16 points total)

### Preparation

``` r
library(VariantAnnotation)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
seqlevelsStyle(txdb) <- "NCBI"
```

### a. Load variant data from a VCF file for `chr8:128,700,000-129,000,000` (2 points)

### b. Reassign the `Seqinfo` for your `vcf` object so that the `Seqinfo` only contains information for chromosomes 1-22 and X. (2 points)

Hint: Use the `keepSeqlevels` function, which can be applied directly to
a `CollapsedVCF` object.

### c. Identify the variants found in promoter regions. (2 points)

Use the `locateVariants` function to find variants overlapping
promoters. Print out the `GRanges` object that is returned by
`locateVariants`.

### d. Combine the fields of the VCF genotype information into a table. (2 points)

If you wish to use `data.table`, refer to *4.1 VCF format* of
*Lecture8\_GenomicData.Rmd* for help.

### e. Retrieve the following information at chr8:128747953. (4 points)

Print out the SNP ID (i.e. “rs ID”), reference base (`REF`), alterate
base (`ALT`), genotype (`GT`), depth (`DP`), allele depth (`ADALL`),
phase set (`PS`).

Hints:

1.  To get the sequence of `DNAString`, use `as.character(x)`.  
2.  To get the sequence of `DNAStringSet`, use
    `as.character(unlist(x))`.  
3.  To expand a list of information for `geno`, use
`unlist(x)`.

### (BONUS) f. Find the phased SNPs for both haplotypes within the region chr8:129127000-129128000. (4 bonus points)
### NOTE: Previously, the region was chr8:128,740,000-128,760,000 but this is a larger region and made the problem much harder. This question is now just for bonus marks.

Print the sequence of phased SNPs for Haplotype 1 and Haplotype 2.
Haplotype 1 is defined as the number to the left of the `"|"` and
Haplotype 2 is defined as the number to the right of the `"|"`. Recall
that a `0` denotes the reference allele and a `1` denotes the alternate
allele.

You may type this out manually in your answer and still the full 4 bonus points.
Bonus question (2 bonus points): provide code to determine the two haplotype
sequences. 
A total of 6 bonus points is available from this question.

**When you are satisfied with your code and answers, use the “Knit”
button in RStudio to create the final set of files that you may then
commit to your private repo and push to GitHub for grading. Please
inspect the files after pushing to GitHub to ensure your code and
figures are rendering as you expect.**
