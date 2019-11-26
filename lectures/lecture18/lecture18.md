## Working with remote data files

### Transferring data

So far, we've used GitHub as a location to store files so we can access them remotely (via `git pull`). This has worked for the small, public files we've been using for class, but this is insufficient for most research projects.

You can use GUI software to transfer data files:
- For transferring data from your computer to a remote storage location, we recommend [Cyberduck](https://sciwiki.fredhutch.org/compdemos/Mountain-CyberDuck/#cyberduck)
- For transferring data among Fred Hutch storage locations, we recommend [Motuz](https://sciwiki.fredhutch.org/scicomputing/store_overview/#data-transfer)

Of course, there are also Unix commands to transfer files, called `scp` (secure copy). More complete instructions can be found [here](https://linuxize.com/post/how-to-use-scp-command-to-securely-transfer-files/), but the most straightforward for our purposes will be to open our shell and transfer a file from our local computer (laptop) to the remote computer (rhino). Navigate to the place on your computer where your chr22 directory is saved (see README of this class for more information) and execute:

  scp chr22.tar USERNAME@rhino:~

Where USERNAME, of course, is your own username, and the filename is accurate as per the file on your computer. The part before the colon (`:`) is the same as the login process to rhino; the part after represents the path to where the file is retained on the cluster. The command above command places the file `setup.sh` in your home directory (`~`) on rhino.

You can copy entire directories by making your command recursive:

  scp -r PATH/TO/DIRECTORY USERNAME@rhino:PATH/TO/RHINO/DIRECTORY

This is especially useful when transferring collections of data files, or entire project directories. You can also use wildcards to specify filenames to transfer.

Reverse the order of the directories to download files from a remote resource:

  scp -r USERNAME@rhino:PATH/TO/RHINO/DIRECTORY PATH/TO/DIRECTORY

If you're downloading files from online, you can use:

  wget URL.to.data.file

This downloads a data file and places it in your current working directory with the same filename as at the end of the URL. If you're trying to download a file in the shell on your local machine and `wget` isn't recognized (this is the case on most Macs), try `curl` instead.

Let's set up a project directory for today's class and try out `wget`:

  mkdir lecture_18
  cd lecture_18
  wget http://genomedata.org/rnaseq-tutorial/practical.tar

### Working with large and/or many data files

If you're working with many data files (especially when transferring among computers), it may be useful to collect them together into a single file. If you receive a file ending in tar (sometimes these are known as tarballs), you can expand them into a regular, accessible directory using:

  tar -xvf DIRECTORY_NAME.tar

Let's try this out with our downloaded data file:

  tar -xvf practical.tar

You can list all the resulting files to see what was in the archive. What kind of files are these?

If the file is also compressed, you can uncompress and untar in a single step by adding one more flag:

  tar -xvzf DIRECTORY_NAME.tar.gz

Let's try this on the file we uploaded via `scp`:

  tar -xvf chr22.tar

Make sure the flag and filename match those on your computer!

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

We can check to make sure our files downloaded appropriately in a few ways. Let's first try to take a peek at the top of one of our data files:

  head hcc1395_tumor_rep1_r1.fastq.gz

Yikes, that doesn't work so well on compressed files! We don't necessarily want to uncompress all our data files, though, so let's try a different tool:

  zcat hcc1395_tumor_rep1_r1.fastq | head

Piping from `zcat` is a great way to inspect sequence files without having to `gunzip` and then `gzip` again.

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

An alternative method for managing environments is `conda`. You have used this on your own computer by installing Anaconda for our Python exercises, and there are methods for working in Anaconda Navigator to [manage environments](https://docs.anaconda.com/anaconda/navigator/tutorials/manage-environments/).

We can also use `conda` on the command line, even in rhino:

  module spider conda

Note that this a method for running Python, too! We can now double check our version:

  conda --version

You can use `conda` to create named environments that will allow you to maintain specific software versions for different projects. For more information on managing `conda` locally, go [here](https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html).

### Setting your path

Let's clear out our loaded modules to make the next section easier to manage:

  module purge

What is in your path? This indicates where your computer will look for executable software:

  echo $PATH

Make a directory to hold the software you want to install (this could be a GitHub repo that contains the scripts you use commonly):

  mkdir ~/my_scripts

Add that location to your path so it is search automatically for executable programs:

  export PATH=$PATH:~/my_scripts

Now you'll see your `my_scripts` folder listed, so that your scripts will automatically be searched when you try to run a command. For more information, see [adding a path to your path](https://opensource.com/article/17/6/set-path-linux).

### Environment variables

The things that start with a `$` in the explanations above represent [environment variables](https://opensource.com/article/19/8/what-are-environment-variables), or information about your session that is used as you execute commands. We can view all environmental variables with:

  env

If we load a new module and then rerun that:

  ml fastqc
  env

...we'll now see FastQC included in our environmental variables.

Note that we can also set our own variables in shell scripts!

## Running software

We talked briefly about viewing and managing processes in our last class. Look [here](https://www.guru99.com/managing-processes-in-linux.html) for more information about these commands:

  ps
  top
  kill

When working on a shared compute cluster, though, we'll need to be a bit more careful about how we run programs (especially large jobs).

### Running jobs interactively

When you log on to the cluster, all commands you run are limited by the compute resources associated with your default login (and by the cluster being shared by many people!). It's good etiquette, and often a necessity, to request your own node for additional tasks:

  grabnode

When prompted, you should type 1 for number of CPUs (since there are lots of us in this class!), hit `Enter` to select the default for memory, and type 1 for number of days. You'll note your command prompt will switch to say `gizmo` instead of `rhino`; you can interact with your directories the same as you would tumorly. See [this documentation](https://sciwiki.fredhutch.org/compdemos/howtoRhino/#guidance-for-use) for more information on when to use `grabnode`). To end your session on `gizmo`:

  exit

This will take you back to your tumor `rhino` login prompt.

### Cluster job submission

In some cases, you may not want to run jobs manually, and will want to use batch computing. In this approach, you specify the parameters for how a job will be run in a script that will then go into a queue. The batch management system on rhino is called [Slurm](http://schedmd.com/). For more information on using slurm on rhino, go [here](https://sciwiki.fredhutch.org/scicomputing/compute_jobs/#using-slurm-on-fred-hutch-systems).

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
- [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) for quality assessment of reads
- [cutadapt](https://cutadapt.readthedocs.io/en/stable/index.html) to remove adapter tags from sequences
- [hisat2](https://ccb.jhu.edu/software/hisat2/manual.shtml) to align/map reads to the reference genome
- [samtools-1.0](http://www.htslib.org/doc/samtools-1.0.html) to manipulate sam files

Let's first do some cleanup on our current data files (make sure you're in `~/lecture_18`):

```
  mkdir data
  mv *fastq.gz data
```

Now we can load the software and assess two of our data files:

```
  ml fastqc
  fastqc hcc1395_tumor_rep1*.fastq.gz
```

What files are output from fastqc? You can view some example output for bad data [here](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/bad_sequence_fastqc.html).

For the following analysis, we're going to proceed with only a single fastq file (to save time here in class).

We have some adapter sequence still included, but we can use cutadapt to remove it:

  ml cutadapt
  cutadapt --adapter=AGATCGGAAGAG --minimum-length=25  -o hcc1395_tumor_rep1_r1_trimmed.fastq.gz hcc1395_tumor_rep1_r1.fastq.gz

We can assume these sequences are clean enough to now continue with analysis.

To map these reads, we need to create a genome index:

  ml hisat2
  hisat2-build chr22.fasta chr22_index

Then we're ready to map the reads:

  hisat2 -p 6 --rg-id=tumor_rep1 --rg SM:tumor --rg LB:hcc1395_tumor_rep1_r1 --rg PL:ILLUMINA --rg PU:H3MYFBBXX.4 -x chr22_index --dta --rna-strandness RF -U hcc1395_tumor_rep1_r1_trimmed.fastq.gz -S ./hcc1395_tumor_rep1.sam

We'll convert the output sam to bam:

  ml samtools
  samtools sort -@ 6 -O sam -T tmp -o hcc1395_tumor_rep1.bam hcc1395_tumor_rep1.sam

From here, we could use the bam file to visualize the genome, assess variants or count reads for differential expression.

## Errata: Making command line work easier

- Configuring your [ssh settings](https://sciwiki.fredhutch.org/scicomputing/access_methods/#advancedoptional-setup-for-making-things-easier-the-ssh-config-file) can save you time
- Questions for discussion:
  - How do we make the analysis above more automated/reproducible?
  - How do you know what programs to use? Which algorithms, and which parameters?
  - What do you do if there's a program you want to use but it's not on the cluster?
  - What if you'll be using a different cluster?
