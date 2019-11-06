# Shell scripting!

We made it!
You have now arrived at the fabled land of shell scripting.
We will slide into shell scripting with the tools we've already built up.


## History revisited

Wouldn't it be nice if we could replay a sequence of commands that we liked?
Here we'll set things up to do exactly that.

Previously we saw that we could view the history conveniently by piping `history` to `less`.
We can see the last set of commands by piping `history` to another command, called `tail`.

Try that!
Look at the manual for `tail` so that you can get the last 25 lines of your history rather than the default.

There's one annoying thing about `history`, which is that it includes the number of the command.
We can avoid that by using [fc](http://pubs.opengroup.org/onlinepubs/9699919799/utilities/fc.html#top) as `fc -ln 1` in place of `history`.
If you get an error with `fc`, use `builtin fc -ln 1`.

Let's noodle around in the shell a bit to generate more history.

* redirect the output of `date` to a file called `date.txt`
* redirect the output of `ls` to a file called `ls.txt`
* echo "shell scripting is hip"
* do your own `echo` command with a phrase of your choice.

Now pipe `fc -ln 1` through tail to get the last 15 commands you have written, and redirect that to a file called `script1.sh`.

Now edit that file with `vi script1.sh`
Delete the lines that you don't want, like `man` calls, etc (an easy way to do that is to go to a line you don't want with the arrow keys, and hit `dd`).
Leave at least one call to `echo` for fun.
Repeat until you are done and then close and save the file `ZZ`.

Now you can execute your shell script with

    bash script1.sh

This runs the commands that you have in `script1.sh`.
You can verify that the script is being run by first deleting the files `date.txt` and `ls.txt` using the `rm` command.
Verify that this removal worked using `ls` then run your script again.

Now it's your turn.
Make another script called `delete.sh` which deletes `date.txt` and `ls.txt`.
Hint, you don't have to start from scratch... you can use your history!


## Making your shell script stand-alone

I hope it's clear that shell scripting is... just a list of commands in a file!
It can get a lot more complicated than that, but it probably shouldn't in most cases.
Once you start needing to do more complex manipulations, it's probably time to reach for a language that has more sophisticated error handling and built-in functionality (such as Python).

You can make your shell script like executing a normal command by

* adding `#!/bin/bash` to the top
* calling `chmod a+x script1.sh`

Now you should be able to run your script with `./script1.sh`.

Next we'll be playing around scripting with samtools using [this bam file](https://console.cloud.google.com/storage/browser/_details/gatk-test-data/wgs_bam/NA12878_20k_b37/NA12878.bam) that's available to you in the GitHub repository.

[samtools](https://www.htslib.org/doc/samtools.html) is a collection of tools for manipulating data in
Sequence Alignment/Map (SAM) format.
This software is installed on rhino,
but needs to be made available for use using `module load`, or `ml` for short:

    ml samtools

Our data has an awkward name, so let's rename it `input.bam` with

    mv wgs_bam_NA12878_20k_b37_NA12878.bam input.bam

This file contains a bunch of sequence reads that align to the human genome.
You first need to prepare this file like so:

    samtools index input.bam

Samtools has a nice command line user interface that has commands within the command.
The format is like

    samtools COMMAND ETC

where `COMMAND` is a command such as `index` and `ETC` are flags, files, etc.
This is now common and is used for git and other lovely software.

You can view reads from the first chromosome with

    samtools view input.bam 1 | head

Try that out.
Rather than using `head`, try piping to something that will allow you to scroll.


## ERRORS!!!

Try running `./script2.sh` via

    ./script2.sh

* What happened?
* Did the script run to completion?
* What is in the final output file?
* Fix the "bug" in this script and re-run.
* How does the output look now?

You see, shell scripts are determined to finish come hell or high water.
That's a really bad thing!
When something goes wrong, such as getting a filename wrong, you want the program to stop, not making bad output.
Now we will learn something very important...

**Always put `set -eu` as the second line of your shell script.**

* Try putting `set -eu` as the second line, just after `#!/bin/bash`.
* Re-introduce the "bug".
* Remove the output files.
* Re-run the script.


## Don't write shell

"Classic" shell is a complex programming language.
The most frequently used shell, called "bash", [is way more complex](https://www.tldp.org/LDP/abs/html/).
However, I do not suggest that you write complex shell scripts.

Rather, I suggest that once you start to need more than just executing a series of commands, you reach for a general-purpose programming language such as Python.
Speaking from personal experience and experience of those in my group, scripts often get increasingly complex to the point where one is doing more complex manipulations that are more suited to a more complete programming language.

However, sometimes we need to execute a command over a series of files, which would typically happen in a loop.
Of course, shell does have loops, but next I'll present an alternative, which is simple and has some good bonuses.


## GNU Parallel

Here we will provide an introduction to a command that provides a convenient and efficient alternative to loops: [GNU Parallel](https://www.gnu.org/software/parallel/).

First, let's load the software we'll need for parallel processing:

    ml parallel

To get our data set up, let's split our input BAM file into a set of smaller files.

    samtools split input.bam

Have a look at these files using your favorite method.

There also exists a command called `file` that tells us information about the file.
Run `file` on `input_0.bam`.

Now, if we want to run file on all of the sam files, we can can use the fact that `file` accepts multiple arguments like so:

    file input_*.bam

If `file` didn't accept multiple files as arguments, we'd have to do something else.

Here we're going to use `parallel` to do that job.
To use it in this case we can do

    ls input_*.bam | parallel file

Here we are piping the file names to the `parallel` command.
We give `parallel` the argument `file`, which tells it which command to run.

Let's try a slightly less trivial example.
There is a program called `gzip` that compresses text files.
We'll play around with that for a bit.

First try

    gzip -k -f input_0.bam

What do the `-k` and `-f` flags do?
Look them up in the `gzip` documentation.

Now try

    ls input_*.bam | parallel gzip -k -f

What did that do?
Use `ls` and `file` to find out.

We now want to start giving flags to `parallel`.
This works like so:

    ls input_*.bam | parallel -v gzip -k -f

This may seem a little confusing at first, but all that's happening is that we are running `parallel` with the `-v` flag, giving `gzip -k -f` as the argument.
What did the `-v` flag do?

OK, now for the fun part: running in parallel!
"Parallel" means to run multiple computations at the same time, which is faster than running them one at a time.
We can tell `parallel` that we want to run (at most) 4 parallel processes at once with the `-j 4` flag:

    ls input_*.bam | parallel -j 4 -v gzip -k -f

Try this version as well as with `-j 1` to run one process at a time.
Which is faster?
Do you notice anything about the execution order when using `-j 4`?

Before we close, I wanted to point out that that the easiest way to run parallel with complex commands is to simply write out a list of those commands to a file, one per line, and just redirect those commands to parallel.
That is, if the commands are `my-commands.sh`, just do

    parallel < my-commands.sh

We have just barely scratched the surface of `parallel`.
It is tremendously powerful!


## More references

* [my notes on shell scripting](https://www.fredhutch.io/articles/2017/02/25/shell-scripting/)
* [GNU Parallel documentation](https://www.gnu.org/software/parallel/).


---

Congratulations!
You are done with the tutorial!
If you want a little more, take a look at the `puritanical.md` file in this same directory.
