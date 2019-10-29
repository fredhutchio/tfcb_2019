# Homework 4: Unix shell

This homework will assess your ability to run commands in the shell, and uses content covered in lecture 09 and 10.

Replace the lines specified in _italics_ with your answers.

The assignment is due on Thursday, November 7 at noon, and includes material from lectures 9 and 12. This assignment will be available in GitHub Classroom on Tuesday, November 5.

## Problem 0

**30 points**

Complete the interactive tutorial.

_Did you hit any points of frustration, and if so, how could we improve the material to avoid that frustration?_


## Problem 1

**25 points**

Learn about the difference between standard out ("stdout") and standard error ("stderr") from [this article](https://www.howtogeek.com/435903/what-are-stdin-stdout-and-stderr-on-linux/) (feel free to read the whole thing, but you can stop before the section "Detecting Redirection Within a Script").
Note that in reading this article, you don't need to come up with a script that will throw an error: we have one at `tfcb_2019/lectures/lecture09/scripting/script2.sh`.

_Write a command here that redirects stdout from `script2.sh` to a file named `stdout.txt` and redirects stderr to a file named `stderr.txt`._


## Problem 3

**10 points**

Use `man` or web search to learn about the `tee` command.

_Modify your previous command to also write stdout to the terminal as well as redirect it to `stdout.txt`_


## Problem 4

**15 points**

Enter the `tfcb_2019/lectures/lecture09/vim` directory.
Have a look at the `anthony.txt` file.
The command

    sed "s/ /\n/g" anthony.txt

will spit out the words of this file one by one.

_Figure out a series of commands, piped together, that will count the number of unique occurrences of various words in `anthony.txt`._

Hint: learn about the `uniq` command and its flags first, but to successfully answer the question you'll need another command to transform your list.


## Problem 5

**20 points**

You might have noticed that the files we're dealing with have "extensions" that describe their file type.
For example, text files are marked with `.txt`, and shell scripts are labeled with `.sh`.

This is a handy convention which is used heavily by a command-line library called "imagemagick" to manipulate images.
Go to the `lecture09/slides/images` directory and try

    convert betty-crocker.jpg betty-crocker.png

which converts `betty-crocker.jpg` (a JPG image) to `betty-crocker.png` (a PNG image).
You can confirm proper conversion using `file`.
Now, your turn:

_Use parallel to convert all of the JPGs in this directory to PNG images._

Big hint: There is a very similar sort of command in the "Compute intensive jobs and substitution" section of the `parallel` man page.

Next:

_Write a script that will take all of the JPGs in the current directory, convert them to PNGs, and then assemble all of the PNGs in the current directory into a file called `montage.png` using the `montage` command. Paste that script here._
