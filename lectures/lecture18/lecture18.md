## Working with remote data files

### Transferring data

So far, we've used GitHub as a location to store files so we can access them remotely (via `git pull`). This has worked for the small, public files we've been using for class, but this is insufficient for most research projects.

You can use GUI software to transfer data files:
- For transferring data from your computer to a remote storage location, we recommend [Cyberduck](https://sciwiki.fredhutch.org/compdemos/Mountain-CyberDuck/#cyberduck)
- For transferring data among Fred Hutch storage locations, we recommend [Motuz](https://sciwiki.fredhutch.org/scicomputing/store_overview/#data-transfer)

Of course, there are also Unix commands to transfer files, called `scp` (secure copy). More complete instructions can be found [here](https://linuxize.com/post/how-to-use-scp-command-to-securely-transfer-files/), but the most straightforward for our purposes will be to open our shell and transfer a file from our local computer (laptop) to the remote computer (rhino). Navigate to `tcfb_2019/lectures/lecture18` on your local computer and then execute:

  scp setup.sh USERNAME@rhino:~

Where USERNAME, of course, is your own username. The part before the colon (`:`) is the same as the login process to rhino; the part after represents the path to where the file is retained on the cluster. The command above command places the file `setup.sh` in your home directory (`~`) on rhino.

You can copy entire directories by making your command recursive:

  scp -r PATH/TO/DIRECTORY USERNAME@rhino:PATH/TO/RHINO/DIRECTORY

This is especially useful when transferring collections of data files, or entire project directories. You can also use wildcards to specify filenames to transfer.

If you're downloading files from online, you can use:

  wget URL.to.data.file.csv

This downloads a data file and places it in your current working directory with the same filename as at the end of the URL. If you're trying to download a file in the shell on your local machine and `wget` isn't recognized (this is the case on most Macs), try `curl` instead.

### Working with large and/or many data files

If you're working with many data files (especially when transferring among computers), it may be useful to collect them together into a single file. If you receive a file ending in tar (sometimes these are known as tarballs), you can expand them into a regular, accessible directory using:

  tar -xvf DIRECTORY_NAME.tar

If the file is also compressed, you can uncompress and untar in a single step by adding one more flag:

  tar -xvzf DIRECTORY_NAME.tar.gz

To create the tarball again, use:

  tar -cvf DIRECTORY_NAME.tar DIRECTORY/

In this command, `DIRECTORY_NAME.tar` is the final output file, which contains the entire contents of `DIRECTORY/`. You can also compress at the same time as tar:

  tar -cvxf DIRECTORY_NAME.tar.gz DIRECTORY/

Please note that `tar` is a notoriously difficult command to remember, resulting in [this classic Unix cartoon](https://xkcd.com/1168/).

You can also use compression algorithms without `tar`. The most common include:

  zip
  gzip
  bgzip

These have their equivalent commands to uncompress as well.

### Confirming data integrity

Compressing and transferring files can sometimes cause problems, especially if the process is interrupted and files are corrupted. You can use [MD5](http://www.digitizationguidelines.gov/term.php?term=md5checksum) (checksums) to check that the entire file is intact. To check a file against another copy, run the following on both files:

  md5sum FILENAME

This creates a unique set of characters that should be identical between the two copies.

## Environments

Your compute environment is the collection of machinery, software, and networks on which you perform tasks. Managing environments is challenging, but there are tools available to help you stay on the right track.

### Finding and loading software

Rhino has hundreds of pieces of software installed, so it uses a module system to make software available for use. To see what software you currently have loaded:

  module list

If you're looking for a specific package, you can search using terms related to the name of the software:

  module spider seq

If you see the package you want, you can load it with:

  module load seqtools

...or the abbreviated version, which is the same as above:

  ml seqtools

Be careful with multiple versions of the same software! You can load specific versions too:

  ml seqtools/4.29

### Conda to control environments

  module spider conda 

### Setting your path

What is in your path?

  echo $PATH

Make a directory to hold the software you want to install (this could be a GitHub repo that contains the scripts you use commonly):

  mkdir ~/my_scripts

Add that location to your path so it is search automatically for executable programs:

  export PATH=$PATH:~/my_scripts

For more information, see [adding a path to your path](https://opensource.com/article/17/6/set-path-linux).

### Environment variables

The things that start with a `$` in the explanations above represent [environment variables](https://opensource.com/article/19/8/what-are-environment-variables), or information about your session that is used as you execute commands. We can view all environmental variables with:

  env

If we load a new module and then rerun that:

  ml fastqc
  env

...we'll now see FastQC included in our environmental variables.

Note that we can also set our own variables in shell scripts!

## Running software

* processes: foreground, background, killing, `top`

### Running jobs interactively

* grabnode

### Cluster job submission

slurm

## Example RNAseq workflow

A standard RNAseq workflow to evaluate genome-wide differential gene expression generally includes the following steps:

1. Biological sampling and library preparation
2. Sequence reads
3. Quality assessment of reads (FastQC)
4. Adapter trimming
5. Map reads to genome
6. Count reads associated with genes
7. Statistical analysis to identify differentially expressed genes

The first two steps happen in the lab, steps 3-6 often occur on the command line in a cluster, and the last step happens in R or Python. The example below works through steps 3-5. For more extended explanations, you can view more complete material [here](https://github.com/hbctraining/Intro-to-rnaseq-hpc-O2/tree/master/lessons) (look in the files with numbered prefixes).

The following code uses these tools:
- fastqc https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

```
ml fastqc
fastqc *.fastq.gz

ml Trimmomatic
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.36.jar PE HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1_paired.fq.gz HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1_unpaired.fq.gz HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2_paired.fq.gz HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2_unpaired.fq.gz ILLUMINACLIP:$EBROOTTRIMMOMATIC/adapters/TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36

ml STAR samtools
chr22.fasta
samtools faidx chr22.fasta
STAR --genomeDir /n/groups/hbctraining/intro_rnaseq_hpc/reference_data_ensembl38/ensembl38_STAR_index/ \
--runThreadN 6 \
--readFilesIn Mov10_oe_1.subset.fq \
--outFileNamePrefix ../results/STAR/Mov10_oe_1_ \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard
```

## Errata: Making command line work easier

* Configuring your [ssh settings](https://sciwiki.fredhutch.org/scicomputing/access_methods/#advancedoptional-setup-for-making-things-easier-the-ssh-config-file) can save you time
