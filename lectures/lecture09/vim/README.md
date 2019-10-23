# Making friends with vi

Vi is a wonderful, powerful, but arcane editor.
It's worth being able to use because even the sparsest linux install will have some variant of vi.
Also, when working with remote machines, it's nice to be able to edit text in a powerful editor directly on that machine.
Finally, sometimes another program (e.g. git) will plop you into vim without you realizing it, so it's nice to know what to do in this situation.
In addition, many programs, such as `less` and Jupyter notebook use vim keybindings for fast navigation.


## The simplest `vi` session

* Invoke `vi` with `vi anthony.txt`
* Type `ZZ` to save and close your file.
* Take a breath.


## Modes

`vi` is a modal editor.
If you have used hotkeys in a web-based email program like gmail you have used modal software before.
For example, if you hit `c` in gmail it will allow you to compose a new email.
When you are typing an email you are in a separate mode-- hitting `c` will just type a `c` in your email.

You can transition from email-writing to hotkey mode by hitting escape (`Esc`) on your keyboard.
`vi` is exactly the same, but the two modes are called "insert mode" and "normal mode."
"Normal mode" is very analogous to hotkey mode but... insanely powerful.

* Invoke vi with `vi anthony.txt`
* Type `i`. Now you are in insert mode. The characters you type are now entered into the document.
* Move around with your arrow keys, and make a few edits, say adding "Hear hear" and "[Applause]".
* When you are done editing, hit `Esc` to exit insert mode and then type `ZZ`
* Take a breath. You are out.


## Colon commands

When you are in normal mode, you also have access to colon commands.

To start you probably only want three of these:

* `:w` saves
* `:q` quits if everything has been saved
* `:q!` quits, even if you have unsaved changes

Now it's your turn!

* Open `anthony.txt` with vim, make a modification, and then save it
* Open `anthony.txt` with vim, make a modification, and then quit without saving


## Moving around quickly

There are lots of ways to move around beyond the arrow keys, but I'm just going to describe three:

* `0` moves to the beginning of the line, and `$` moves to the end of the line
* `b` moves back one word, and `w` moves forward one word
* `{` moves back one paragraph, and `}` moves forward one paragraph

You can prefix these commands with numbers to move faster, e.g. `3w` moves you forward three words.

`Ctrl-o` jumps you back to previously visited locations, while `Ctrl-i` jumps you forward.
You can think of these as undo and redo for movement.

Take some time to move around `anthony.txt`, making some text modifications as you go.


## Search

To search, use `/` and type what you want to find.
When you are done typing the thing you want to find, hit return.
If you want to find the next item, hit `n`.
You can go back using `Ctrl-o` as described before.


## Cutting and pasting

There is a simple way to cut and paste using vim which is exactly analogous to a word processor: highlight a block of text, then copy or cut, then paste.

* Move to where you want to start your highlight
* Press `v` (this places you in "visual mode")
* Move to the end of your highlighted region
* Press `d` to cut, `y` to copy
* Move to where you want to paste
* Press `p` to paste

You can also cut and paste using `d` and `y` together with a motion key (e.g. `dw` cuts a word, and `d2w` cuts two); `dd` and `yy` cut or copy an entire line, respectively.

Try doing some cutting and pasting.


## Undo/redo

If you ever mess anything up (which is easy to do in command mode), `u` is undo and `Ctrl-r` is redo (from command mode).


## Everyone does this once

To save, do *not* use `Ctrl-s`.
This will only lock your terminal, which you then must rescue by using `Ctrl-q`.
As described above, use escape to leave insert mode then call `:w`, or `:wa` to save all files you might have open.

Relatedly, if you hit `Ctrl-z` in normal mode, it will pause the program and return you to the shell.
To get back to your `vi` session, type `fg`.


## Vim resources

Here we'll be a little specific and say that the software we've been calling `vi` is actually a modern variant called [vim](https://www.vim.org/), which is written by a Dutch fellow who uses the software to raise money for children in Uganda.

* a [cheat-sheet](http://i.imgur.com/YLInLlY.png)
* a [wallpaper](https://github.com/LevelbossMike/vim_shortcut_wallpaper)
* [Ben Crowder's vim tips](http://bencrowder.net/files/vim-fu/)
* The `vimtutor` command, available wherever you find `vim`
* [An online vim tutorial](http://www.openvim.com/)

Note that these guides insist that you can't use the arrow keys.
You certainly can, though it's not considered hip (because vim is all about efficiency, and moving your hands from home position to the arrow keys is not efficient.)


---

Congratulations!
You now know a little `vi` which we'll use for the rest of this course.

Head over to the `lecture09/history` directory for the next set of instructions.
