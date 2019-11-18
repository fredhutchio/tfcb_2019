# Lecture 17: Introduction to remote computing

**Erick Matsen ([@ematsen](https://twitter.com/ematsen), [matsen.group](http://matsen.group))**

Now that you have some experience with shell, we'll dive a little more deeply into what remote and cloud computing are, including the modern container-based approach to high-performance computing.
However, we will start by going over some of the questions from last time.

## Learning objectives

After this course, you should be able to:

- Be able to run shell commands from Python
- Understand local versus remote execution
- Understand the basics of what cloud computing is and what it offers
- Understand what a container is and how

## Class materials

This directory has the class materials for the course.
I have added a `quickref.md` quick reference file with basic commands we've used so far.

## Reminders


## Interactive work

### A little about wildcards/globbing

In the scripting component from last time, we used wildcards like so:

    ls input_*.bam | parallel file

The `*` just represents any set of characters.
Thus this `ls` call just lists anything that starts with `input_` and ends with `.bam`.

Now it's your turn.
You can combine wildcards, using them for example for both directory and file names.
Try using a wildcard to list all of the `.md` files in `../lecture09/directories/purchase`.

There is also a Python `glob` module which implements wildcards.
Here's an example you can copy and paste into the shell:

    python -c "import glob; print(glob.glob('ps*'))"

Here we give the `glob.glob` command a string to expand, which in this case is like using `ls ps*`.
`glob.glob` returns a Python list.

[There are other symbols you can use if you want more precise control over your wildcard matches](http://www.linfo.org/wildcard.html) but `*` is a fine start.


### Hidden directories and files

The shell has a special way of hiding files that you don't want to see: any file or directory starting with a period (`.`) is hidden when you `ls`.
It's also not matched with `*`.
Try doing

    echo ~/*

which compactly lists all of the non-hidden files in your home directory (which is `~`).
Then try

    echo ~/.*

which shows all of the hidden files.

You can also use

    ls -a ~

which will list all of the files in `~`.


### Removing files

Removing files happens with

    rm FILENAMES

where `FILENAMES` can be 1 or more names.
Removing directories happens with

    rmdir DIRECTORYNAMES

where `DIRECTORYNAMES` is 1 or more directory names.
However, these directories need to be empty.
It is more common to use the less safe

    rm -r DIRECTORYNAMES

which recursively removes `DIRECTORYNAMES` and all of their contents.

Unlike "putting something in the trash" on a desktop OS, once linux removes a file, it's gone!

However, there's a safety net on the Fred Hutch linux system, which are the "snapshots."
Every hour the Fred Hutch servers back up your whole home directory into a snapshot folder hosted in `~/.snapshot`.
Try

    ls ~/.snapshot

to see all of these files.

Now it's your turn to pretend you made a mistake.
I will assume that you are in `tfcb_2019/lectures/lecture17`.
Try removing some non-essential file, say doing

    rm ../lecture09/vader.txt

and then restoring it from the snapshot directory.


### Making directories

Making directories is done with

    mkdir DIRECTORY

which makes a directory called `DIRECTORY`.


### Copying and moving/renaming files

Copying files uses the `cp` command, like so

    cp SOURCE DESTINATION

where `SOURCE` is an extant file, and `DESTINATION` is a place the file can go.
If `DESTINATION` is a file name, it will copy that file to have that new name.
If `DESTINATION` is a directory, it will copy the file with the same name into that directory.

You can also use the form

    cp SOURCES DIRECTORY

where `SOURCES` is a collection of files and `DIRECTORY` is a directory.

Moving/renaming files is just the same, except that you use the `mv` command rather than the `cp` command.

Now it's your turn.

* Copy this `README.md` file into `..`
* Rename it `../README.copy.md`
* Make a directory called `copies`
* Copy `../README.copy.md` and `../lecture09/vader.txt` into `copies`
* Remove the `copies` directory


### Trying out mixing shell and Python

The `file` command, when run on `../lecture09/slides/images/betty-crocker.jpg`, returns

    ../lecture09/slides/images/betty-crocker.jpg: JPEG image data, JFIF standard 1.01, resolution (DPI), density 72x72, segment length 16, baseline, precision 8, 460x600, components 3

which as you can see, is the file name, then a colon, then a description of the file.

Make a Python script that runs `file` on each file in `../lecture09/`, extracts just the descriptions, and sorts them.

You can use the `ps-count.py` script in this directory as a template (in particular the invocation of `subprocess.check_output`).
