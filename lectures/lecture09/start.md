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

Many commands have built-in help using `-h` or `--help`.
Try `less --help` to find the commands that allow you to move forward and backward by page, and use them to view `sequence.gb`.

For more in depth help see the command `man`, which is short for "manual."

OK, so now try `man man` to learn about the command `man`.

1. first just enter `man` and then quit (also with the letter `q`)
2. then re-enter `man` and use the online help to figure out how to forward and backward one window at a time
   (note that these same commands work for the `less` command)
3. now look through the contents of the `man` manual to find a flag that will show you the location of the file containing the documentation
4. save this information somewhere (on a piece of paper?)
5. execute `cd directories` then `less README.md` to work through the next part of the adventure
