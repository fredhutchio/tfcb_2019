## Working with remote data files

### Transferring data

So far, we've used GitHub as a location to store files so we can access them remotely (via `git pull`). This has worked for the small, public files we've been using for class, but this is insufficient for most research projects.

You can use GUI software to transfer data files:
- For transferring data from your computer to a remote storage location, we recommend [Cyberduck](https://sciwiki.fredhutch.org/compdemos/Mountain-CyberDuck/#cyberduck)
- For transferring data among Fred Hutch storage locations, we recommend [Motuz](https://sciwiki.fredhutch.org/scicomputing/store_overview/#data-transfer)

Of course, there are also Unix commands to transfer files, called `scp` (secure copy). More complete instructions can be found [here](https://linuxize.com/post/how-to-use-scp-command-to-securely-transfer-files/), but the most straightforward for our purposes will be to open our shell and transfer a file from our local computer (laptop) to the remote computer (rhino). Navigate to  `tcfb_2019/lectures/lecture18` on your local computer and then execute:

  scp setup.sh USERNAME@rhino:~

This places the file `setup.sh` in your home directory (`~`) on rhino.

You can copy entire directories by making your command recursive:

  scp -r PATH/TO/DIRECTORY USERNAME@rhino:PATH/TO/RHINO/DIRECTORY

  wget
  curl

### Working with large and/or many data files

* uncompressing zip files from the command line

  tar

  zip
  gunzip

## Environments

### Environment variables

What is an environment?

* environment variables; PATH

### Finding and loading software

  module list

  module spider

  module load
  ml

## Running software

* processes: foreground, background, killing, `top`

### Running jobs interactively

* grabnode

Example bioinformatics workflow
- fastqc https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
- Trimmomatic http://www.usadellab.org/cms/?page=trimmomatic

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

### Cluster job submission

slurm

## Making command line work easier

* ssh access keys
* ssh config https://sciwiki.fredhutch.org/scicomputing/access_methods/#advancedoptional-setup-for-making-things-easier-the-ssh-config-file
