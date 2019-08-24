# Directories

Congratulations, you just changed directory using `cd`, which is short for "change directory".
By executing `cd directories` you moved your current location inside a directory called `directories`.
You can think about this being like double clicking on a folder in a graphical user interface.

Think back to the previous lesson, and execute the command that allows you to see the currently available set of files from this directory.

### Paths

It takes a little work to keep track of where you are in the shell.
Try executing the command

    pwd

to show where you are.
The `pwd` command is short for "print working directory."

On my computer the `pwd` command outputs

    /home/ematsen/writing/tfcb_2019/lectures/lecture09/directories

which shows the list of directories I'm in.
Such a list of directories (perhaps ending with a file name) is called an _absolute path_.
Each successive directory is separated from the previous one with a `/`, so I'm in a directory called `directories`, which is contained in a directory called `lecture09`, which is contained in a directory called `lectures`, and so on.
As I described before the beginning ("root") directory is called `/`.

Remember how before we learned how to give an argument to the directory-listing command?
Now give that command a part of your path, so see the list of files at various places in your directory hierarchy.

### Tree

Of course, directories can be nested and so it's not always easy to get an overall view of what is where.
For this there is the `tree` command, which will show a nested directory structure as a tree:

    .
    ├── purchase
    │   ├── egg
    │   │   └── yolk.md
    │   └── greenland
    │       └── ice.md
    └── README.md

    3 directories, 3 files

Try executing `tree` now and see if you get the same result.

The `tree` command can also take an argument, which is the directory for which it should draw the tree.

* without changing directories, show the directory tree of the `purchase` directory only
* look in the `tree` documentation to list only directories and try that out here

### Moving around more freely

Try `cd`-ing into the `purchase` directory.
How do we get back?

First try the slow way.
Use `pwd` to find your current directory as an "absolute path" starting from the root directory `/`.
Then copy and paste the that whole path, except for the `purchase` at the end, and use that as an argument for `cd`.
Use `pwd` again to check your location.

There is an easier way!
The `..` symbol means "one directory towards the root," so `cd ..` will take you to the directory that contains your current directory.
Try that out now.

We can also `cd` through multiple directories at once by giving an argument like `cd first/second/third` (replacing `first`, `second` and `third` with actual names of directories.
Try jumping into the `greenland` directory now.
