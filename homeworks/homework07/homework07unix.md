# Homework 7: Bioinformatics with remote computing

Log in to rhino and use data found at http://genomedata.org/rnaseq-tutorial/HBR_UHR_ERCC_ds_5pc.tar to answer the following questions.

For each question, replace the placeholder text (for `code` and _in italics_) with your answers.

## Problem 0

**5 points**

Transfer the data to your home directory on rhino and unarchive/expand the data files. How many files are included in the original file, and what is the file format?

```
Replace this text with your code
```

_Report the number of files and file format here_

## Problem 1

**5 points**

What is the sequence identifier of the 100th sequence in the file? Hint: the line containing the identifier for each sequence starts with `@`, and there are four lines for each sequence.

```
Replace this text with your code
```

_Report the sequence identifier for the 100th sequence here_

## Problem 2

**10 points**

The dataset for the python portion of this homework is relatively small. If you were given a much larger dataset and needed to begin testing similar code to run on a remote cluster, what would you need to do to ensure the software was available to run on rhino? Note: assume you are running a python script, not an ipython notebook.  

```
Replace this text with your code
```

_Explain the code above here_

## Problem 3

**10 points**

In lecture 18, we ran bioinformatics command-line software on a single data file from our dataset. In a real research project, we may have dozens (or even hundreds) of data files.

Describe two ways you could run `fastqc` on multiple data files by executing code only once (e.g., without having to enter `fastqc FILENAME` for every file or sample). What are the advantages and disadvantages of each method? Also relate how each of those approaches would be more complicated for the code we used to execute `cutadapt`.

_Report your answer here_
