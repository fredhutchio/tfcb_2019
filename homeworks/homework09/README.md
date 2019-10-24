# Homework 9

This homework will assess your ability to run commands in the shell, and uses content covered in lecture 09.

Please make a copy of this file, and replace the lines specified in _italics_ with your answers then submit your homework to GitHub Classroom via a new repository containing that file.


## Problem 0

**40 points**

Complete the interactive tutorial through the section on redirection.

_Ask a question here._


## Problem 1

**25 points**

Learn about the difference between standard out ("stdout") and standard error ("stderr") from [this article](https://www.howtogeek.com/435903/what-are-stdin-stdout-and-stderr-on-linux/) (feel free to read the whole thing, but you can stop before the section "Detecting Redirection Within a Script").
Note that in reading this article, you don't need to come up with a script that will throw an error: we have one at `tfcb_2019/lectures/lecture09/scripting/script2.sh`.

_Write a command here that redirects stdout from `script2.sh` to a file named `stdout.txt` and redirects stderr to a file named `stderr.txt`._


## Problem 3

**10 points**

Use `man` or web search to learn about the `tee` command.

_Modify your previous command to also write stdout to the terminal as well as redirect it to `stderr.txt`_


## Problem 4

**15 points**

Enter the `tfcb_2019/lectures/lecture09/vim` directory.
Have a look at the `anthony.txt` file.
The command

    sed "s/ /\n/g" anthony.txt

will spit out the words of this file one by one.

_Figure out a series of commands, piped together, that will count the number of unique occurrences of various words in `anthony.txt`._

Hint: learn about the `uniq` command and its flags first, but to successfully answer the question you'll need another command.


# Problem 5

**10 points**

The `grep` command can act as a filter that you can use as part of your series of piped commands.
For example, `grep "f"` will only select words that contain "f".
Read `grep`'s documentation and find flags so you can...

_Modify your answer to the previous question to obtain the number of unique occurrences of words that do NOT contain uppercase or lower-case f in a `anthony.txt`._
