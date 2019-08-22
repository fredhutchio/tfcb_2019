# Lecture 9

     _______     _______.  ______     ___      .______    _______
    |   ____|   /       | /      |   /   \     |   _  \  |   ____|
    |  |__     |   (----`|  ,----'  /  ^  \    |  |_)  | |  |__
    |   __|     \   \    |  |      /  /_\  \   |   ___/  |   __|
    |  |____.----)   |   |  `----./  _____  \  |  |      |  |____
    |_______|_______/     \______/__/     \__\ | _|      |_______|

    .___________. __    __   _______
    |           ||  |  |  | |   ____|
    `---|  |----`|  |__|  | |  |__
        |  |     |   __   | |   __|
        |  |     |  |  |  | |  |____
        |__|     |__|  |__| |_______|

         _______. __    __   _______  __       __
        /       ||  |  |  | |   ____||  |     |  |
       |   (----`|  |__|  | |  |__   |  |     |  |
        \   \    |   __   | |   __|  |  |     |  |
    .----)   |   |  |  |  | |  |____ |  `----.|  `----.
    |_______/    |__|  |__| |_______||_______||_______|


## Starting out

To get started, simply type

    ls

and hit return.

This shows the collection of files that are in your current directory (a.k.a. folder).
The name `ls` is an abbreviation for "list".

This is the simplest type of shell command invocation.


## Anatomy of a shell command

More generally, shell commands look like

    command flags arguments

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

    ls -w 15 /

in which we give the `-w` flag the value of 15.
The `-w` flag tells `ls` that the user wants to wrap output to some number of letters wide and no wider, in this case 15 letters wide.

Try varying the argument to `-w` from very small to very large.
Are negative numbers allowed?

All flags start with a dash , like `-1` and `-w` did.
They can be longer than one letter, such as `-name`.

Some flags start with two dashes, like `--verbose`.
Generally two-dash flags are for long names and single-dash flags are for short names, but this convention isn't always consistently applied.

For example, try using `--width` in place of `-w` in the command above.
Does `-width` work?

---

Congratulations, you know know the basics of command invocation in the shell.

To start your adventure, type `cat start.txt` from inside the `lecture09` directory in the TFCB repository.
