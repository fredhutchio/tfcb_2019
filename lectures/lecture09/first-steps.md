# First steps

To get started, simply type

    ls

and hit return.

This shows the collection of files that are in your current directory (a.k.a. folder).

The name `ls` is an abbreviation for "list".

This is the simplest type of shell command invocation.

Next just mash the keyboard making a random command, like this:

    asdfasfkjaeasdkjaadsfkljaeflkwae

What happens?


## Anatomy of a shell command

More generally, shell commands look like

`command` `flags` `arguments`

where `flags` and `arguments` are optional.
We separate these components with spaces (extra spaces are fine).


## Arguments

Try this command with an argument:

    ls /

Here we are passing `/` to the `ls` command as an argument to tell it which directory should be listed.
Specifically, we are telling it to list the very "bottom" directory of the file system, called root.
We will discuss directory structure later.

Now try

    ls asdfasdf

which should give an error unless you happen to have a file named `asdfasdf` in your current directory.
Errors are OK.
They don't cause any harm, and are a natural part of working in the shell.


## Flags

Try this command with a flag:

    ls -1 /

Here `-1` is a flag which specifies that the listing should happen in a single column.
Flags generally go before the arguments, but some commands will allow them in any order.

Does `ls` allow flags after arguments?
Try it out!

Flags modify some form of the command behavior, and can range from simple to rather complex.

Flags can have arguments, like

    ls -w 50 /

in which we give the `-w` flag the value of 50.
The `-w` flag tells `ls` that the user wants to wrap output to some number of letters wide and no wider, in this case 50 letters wide.

Try varying the argument to `-w` from very small to very large.
Are negative numbers allowed?

All flags start with a dash , like `-1` and `-w` did.
They can be longer than one letter, such as `-name`.

Some flags start with two dashes, like `--verbose`.
Generally two-dash flags are for long names and single-dash flags are for short names, but this convention isn't always consistently applied.

For example, try using `--width` in place of `-w` in the command above.
Does `-width` work?

Try the `-a` flag.
Now try the `-a` and `-1` flags together as `ls -a -1`.


## Looking at files

Try

    cat vader.txt

This is a handy way to look at small files.
Now try

    cat sequence.gb

This is not a very convenient way to view a long file, because it spits out the whole file and you will no longer be able to see the top.

A better way is to use `less`.
Before we try that, note:
**you can navigate in `less` using the arrow keys, and exit using the letter `q`**.
Try

    less sequence.gb

and have a look around this file.
Exit that `less` session using `q`.

Many commands have built-in help using `-h` or `--help`.

Try

    less --help

to find the commands that allow you to move forward and backward by page, and use them to view `sequence.gb`.


## `man`: the built-in manual

For more in depth help see the command `man`, which is short for "manual."
The argument to `man` is the command you want to see documentation for.

Try

    man man

to learn about the command `man`.

0. Look at the overall structure of the document. What does the top describe? What about the rest?
1. now look through the contents of the `man` manual to find a flag that will show you the location of the file containing the documentation
2. save this information somewhere (on a piece of paper?)


## On your own

Use any method you like to find documentation about the `wc` command, then use it to count the number of lines in `sequence.gb`.


---

Congratulations!
You now know how to navigate around the file system.

For the next step, look at the README file in the `directories` directory to work through the next part of the adventure.
