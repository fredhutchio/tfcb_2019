When you ran `cat start.md` it spat the contents of this file out to the shell.
This is a handy way to look at small files.
Now try

    cat sequence.gb

This is not a very convenient way to view a long file, because it spits out the whole file and you will no longer be able to see the top.

A better way is to use `less`.
Before we try that, note:
**you can navigate using the arrow keys, and exit using the letter `q`**.
Try

    less sequence.gb

and have a look around this file.

We are now going to learn about the command `man`, which is short for "manual."

OK, so now try `man man` to learn about the command `man`.

1. first just enter `man` and then quit
2. then re-enter `man` and use the online help to figure out how to forward and backward one window at a time
   (note that these same commands work for the `less` command)
3. now look through the contents of the `man` manual to find a flag that will show you the location of the file containing the documentation
4. save this information somewhere (on a piece of paper?)
5. execute `cd directories` then `cat README.md` to get to the next part of the adventure
