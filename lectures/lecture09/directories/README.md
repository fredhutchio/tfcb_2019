# Directories

Congratulations, you just changed directory for the first time.
By executing `cd directories` you moved your current location inside a directory called `directories`.
You can think about this being like double clicking on a folder in a graphical user interface.

Think back to the previous lesson, and execute the command that allows you to see the currently available set of files from this directory.

It takes a little work to keep track of where you are in the shell.
Try executing the command

    pwd

to show where you are.
The `pwd` command is short for "print working directory."

On my computer the `pwd` command outputs

    /home/ematsen/writing/tfcb_2019/lectures/lecture09/directories

which shows the list of directories I'm in.
Such a list of directories (perhaps ending with a file name) is called a _path_.
Each successive directory is separated from the previous one with a `/`, so I'm in a directory called `directories`, which is contained in a directory called `lecture09`, which is contained in a directory called `lectures`, and so on.
As I described before the beginning directory is called `/`.

Remember how before we learned how to give an argument to the directory-listing command?
Now give that command a part of your path, so see the list of files at various places in your directory hierarchy.

introduce tree
..

    .
    ├── purchase
    │   ├── egg
    │   │   └── yolk.md
    │   └── greenland
    │       └── ice.md
    └── README.md

    3 directories, 3 files
