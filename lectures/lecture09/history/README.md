# History

Your "history" is the list of commands you have entered.

## Navigating your history

The easiest way to get back in your history is just to hit the up arrow.
Repeatedly hitting up arrow will take you back in time, while down arrow takes you forward in time.
When we see a command that we like, we can hit `Return` to execute that command, or use the left and right arrow keys to move to a place where we can edit it.

Here are some commands to help you browse that history:

* `history`: gives the full history; likely to be too much too be useful
* `history | less`: history made more navigable

Now your turn!
Take a look through your history.

You can also search through your history.
Hitting `Ctrl-r` brings up reverse interactive search.
Here you could type `cd`, which brings up the most recent command containing the string `cd`:

```
$ cd history
bck-i-search: cd_
```

You can cycle through earlier commands by hitting `Ctrl-r` again.
If you want to cancel your reverse search, use `Ctrl-g`.

Try looking back through your history to find all the commands that used `tree`.


## Command line editing

Maybe you just put together a big command, but realized that there's some problem at the very beginning.
Do you use the arrow keys to go one letter by letter to the start of the line?
No!
You use the following handy commands:

* `Ctrl-a`: move to the start of the line
* `Ctrl-e`: move to the end of the line
* `Meta-f`: move forward one word
* `Meta-b`: move backward one word
* `Ctrl-k`: delete to the end of the line

What is `Meta`?
If you are on a PC keyboard, it's the `Alt` key.
If you are on a Mac, it may be `Option` already, or [a little configuration may be required](http://osxdaily.com/2013/02/01/use-option-as-meta-key-in-mac-os-x-terminal/).

On the other hand, if you are using Ubuntu (e.g. you are using the Fred Hutch servers) then you should just be able to use `Ctrl` with the left and right arrows to move by word.

Another option for you `vim` lovers out there is the `fc` command, for "fix command."
`fc` puts your last command into a little `vim` session, allowing you to modify it as you like.
When you save and quit it executes the modified command.
If you decide that you don't want to do anything, you can just remove all text from the file.


## For OS X users who want a long history.

It appears that OS X truncates your history at a measly 500 lines.
Phooey on that!
Put this in your `~/.bashrc` file to get an unlimited history:

```
export HISTFILESIZE=
export HISTSIZE=
```

More about your `.bashrc` below.
Note that [this isn't a perfect solution and a better one exists](http://superuser.com/a/664061), but it's good enough.


---

Congratulations, now you know something about command history and a little more about command line editing.
We're now ready to do a little shell scripting!
Next navigate to `lecture09/scripting`.
