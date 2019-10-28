# Lecture 8: Genomic data in R

This lecture will unite the last lecture's content on genomic analysis with our previous coding in R. The packages we'll use this week are from [Bioconductor](http://bioconductor.org), a collection of software specifically designed for genomic analysis in R.

## Learning objectives

- Use Bioconductor packages to work with genomic data in R
- Load, inspect, and query genomic data (BED/SEG, BAM, VCF files)
- Identify and annotate genomic variants

## Class materials

- If you have not done so already, update your local copy of the class repository from GitHub. You should have a directory (`lecture08`) containing the following three RMarkdown tutorials:
  - [Lecture8_GenomicData.Rmd](Lecture8_GenomicData.Rmd): storing genomic data as objects, assessing genomic ranges, importing and assessing variant (vcf) data
  - [Lecture8_Annotations.Rmd](Lecture8_Annotations.Rmd): load and apply gene annotations to genomic data
  - [Lecture8_Rsamtools.Rmd](Lecture8_Rsamtools.Rmd): Compute “pile-up” statistics at genomic loci to identify genomic variants
- Please download all data files found in [this folder](https://www.dropbox.com/sh/zoitjnobgp7l7c2/AABBIpTQcNA4lWYOFnV5dlMKa?dl=0) and add them to your `lecture08` directory. The files should have the following filenames:
  - `BRCA.genome_wide_snp_6_broad_Level_3_scna.seg`
  - `BRCA_IDC_cfDNA.bam`
  - `BRCA_IDC_cfDNA.bam.bai`
  - `GIAB_highconf_v.3.3.2.vcf.gz` (if this file was automatically uncompressed on your computer, resulting in a file named `GIAB_highconf_v.3.3.2.vcf`, look in your Trash folder to find the original file ending in `gz`)
  - `GIAB_highconf_v.3.3.2.vcf.gz.tbi`
- You should run [this script in RStudio](../../software/genomic_data.R) to ensure all Bioconductor packages are installed.
- Your homework will also require use of the [Integrative Genomics Viewer (IGV)](https://software.broadinstitute.org/software/igv/)

## Reminders

- Homework 2 was due at noon on Tuesday, October 22 in GitHub Classroom.
- Homework 3 (genomic data in R) is available through GitHub Classroom and is due Tuesday, October 28 at noon. You should receive an email containing an invitation to create your repository using GitHub Classroom. Contact Kate (khertwec at fredhutch.org) with any questions or concerns.
- The next class session will cover Unix command line. Please ensure you have the [software for accessing the Unix command line](https://github.com/fredhutchio/tfcb_2019/tree/master/software#unix-command-line) installed. We will be using compute clusters available through Fred Hutch, which will require use of your HutchNetID. If you have not yet received a notification of HutchNetID creation, please contact Kate.
