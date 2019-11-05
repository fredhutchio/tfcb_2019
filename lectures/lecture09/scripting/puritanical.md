# Puritanical shell scripting

This post and talk are inspired by Ryan Tomayko's [talk](http://confreaks.tv/videos/gogaruco2010-the-shell-hater-s-handbook) *The Shell Hater's Handbook*.

The "puritanical" philosophy of shell scripting is that shell is *not* a general-purpose programming language.
Shell is a programming language that is designed for assembling and running other commands.
Use it for this task, and this task only.

Puritans know that if you need something more than composition of commands, it's time to reach for a more suitable programming language such as Python.
Python will be able to handle your needs as they expand and become more complex.

We will be describing POSIX shell, which a standardized set of commands that are available in any modern standards-compliant shell such as `bash`, `zsh`, and `dash`.
Note that bash adds lots of features on top of the POSIX standard, such as arrays, but again if you are feeling the need for these features my opinion is that a more general-purpose programming language is a better choice.
Nevertheless, shell scripts are very handy and a great tool.

The best reference manual for POSIX shell is [this handbook](http://shellhaters.org/), which indexes directly into the POSIX documentation.
For a more introductory perspective, try [Classic Shell Scripting](http://shop.oreilly.com/product/9780596005955.do) by Robbins and Beebe (you *may* be able to find a free PDF of this book using your search engine).

* [Your first shell script](#your-first-shell-script)
* [Using predefined variables](#using-predefined-variables)
* [Quoting](#quoting)
* [Command substitution](#command-substitution)
* [Defining variables](#defining-variables)
* [Return codes](#return-codes)
* [Control flow](#control-flow)
* [The tiniest bit of `sed`](#the-tiniest-bit-of-sed)
* [Functions](#functions)
* [Command line arguments](#command-line-arguments)


## Your first shell script

Shell scripting is to first approximation just writing what you would write in the interactive shell in a file and then running it.
`echo` is a command that simply echoes its arguments back to you.

```
#! /bin/sh

set -eu

echo "Twinkle twinkle little *"
```

Save this in `script.sh` and make it executable using `chmod +x script.sh`.
Now we can execute it:

```
$ ./script.sh
Twinkle twinkle little *
$
```

As you can see, to execute a command in a shell script, you just put it on a line and it gets called.
[ *Important note:* we're using a `$` at the beginning of a line to denote your shell prompt. You shouldn't type this in! I'll only use it when showing output. ]

Now, what is this `set -eu` bit?
Although the script will run without it, I strongly suggest that every script you write contain this line.
Adding it will make your shell script fail and print an error message when either a command fails or a variable substitution did not work.
That's what you would expect, right?

You might also consider `set -eux`, which prints every executed line.


## Using predefined variables

One of the key aspects of shell is variable substitution: substituting in the value of a variable.
To demonstrate, we'll use some pre-defined variables.
We can do this directly in the interactive shell.

```
$ echo "My name is $USER and my editor is $EDITOR."
My name is matsen and my editor is vim
$
```

There are lots of pre-defined variables, and you can see them using the `env` command.


## [Quoting](http://pubs.opengroup.org/onlinepubs/9699919799/utilities/V3_chap02.html#tag_18_02)

If you don't know about quoting, there is a key point here: the quotes make the argument of `echo` a single item.
If instead we tried

```
echo Twinkle twinkle little *
```

this echoes `Twinkle twinkle little` with all of the other files I have in my directory, because `*` is a pattern that [matches all characters](http://pubs.opengroup.org/onlinepubs/9699919799/utilities/V3_chap02.html#tag_18_13).

You might have noticed that the double quotes with variables above didn't actually perfectly preserve the string we put into echo.
Rather, the value of the variables were substituted in.
Sometimes this is the effect we want, but sometimes we literally want a string with a dollar sign.

For this there are single quotes:

```
$ echo 'My name is $USER and my editor is $EDITOR'
My name is $USER and my editor is $EDITOR
$
```

Note that we could have achieved the same effect using backslashes, which remove the special meaning of characters.

```
$ echo "My name is \$USER and my editor is \$EDITOR"
My name is $USER and my editor is $EDITOR
$
```


## Command substitution

We'd like to be able to use the result of one command as the input for another.
We can do this using backticks (i.e. wrapping your command in \`'s), but I prefer the more explicit `$( )` notation, which also nests properly.

For example, say we have a nested directory structure like so:
```
$ tree
.
├── buckle
│   └── shoe
└── one
    └── two
```

When we `ls`, we get the first level:

```
$ ls
buckle  one
```

When we call `ls $(ls)`, we are calling `ls` on the result of calling `ls`, which as we saw above is `ls buckle shoe`:

```
$ ls $(ls)
buckle:
shoe

one:
two
```

_Exercise:_ Learn about the `sleep` command, and write a command that sleeps as many seconds as there are lines in the `script2.sh` file.
Hint: if you directly call `wc` on a file it will list the file that it is called on as part of the output, however if you pipe `script2.sh` to `wc` it will not.


## Defining variables

You know now that variables are there to store strings for later substitution.
You define them like so:

```
color="chicken"
```

Note that there can be no spaces around the equals sign when assigning variables.

It's common to use all caps for shell variables:

```
COLOR="chicken"
```

which means that you won't confuse them with commands, though you don't have to follow that convention.

Now we can use the variable:

```
$ COLOR="chicken"; echo "Roses are red, violets are $COLOR"
Roses are red, violets are chicken
```
(I snuck in a semicolon here, which allows us to execute several commands on a single line.)
It's common to set variables to the output of commands.
For example, try this:

```
DATE=$(date); sleep 3; echo "3 seconds ago it was $DATE"
```

So far we've been using variables in situations where the termination of the variable is clear, such as at the end of the string.
This need not always be the case.
For example, we might want to substitute in the value of `COLOR` as so:

```
COLOR="blue"; echo "Let's paint the house $COLORish"
```

Instead of getting the output "Let's paint the house blueish", the shell interprets `COLORish` as its own variable, which isn't defined.

The right syntax is `${VARIABLE}otherstuff`:

```
$ COLOR="blue"; echo "Let's paint the house ${COLOR}ish"
Let's paint the house blueish
$
```

_Exercise:_ echo a string that describes your status, including these variables: `USER`, `PWD`, and `SHELL`.


## Return codes

Return codes are how programs communicate if they succeeded or how they failed.
They also form the basis of the conditional execution system in shell.
In principle they go from 0 to 255, but you only need to know that 0 is success and anything else is a failure.
This might sound backwards (and is indeed opposite of languages such as Python), but in shell there are lots more ways to fail than there are to succeed.
We have heard of these concepts above in the discussion of `set -eu`.

The return code is held in the `$?` variable.

```
echo "the most recent return code is $?"
this_command_doesnt_exist
echo "the most recent return code is $?"
```

These concepts are used in AND lists and OR lists.

* `cmd && other_cmd`: `other_cmd` will only be run if `cmd` succeeds
* `cmd || other_cmd`: `other_cmd` will only be run if `cmd` fails

Try this:

```
this_command_doesnt_exist || echo "the previous command failed"
```

More return code action:

```
$ ls gefilte_fish && echo "we have fish"
ls: cannot access 'gefilte_fish': No such file or directory
$  ls gefilte_fish || echo "we have no fish"
ls: cannot access 'gefilte_fish': No such file or directory
we have no fish
```

We only got the echo message in the second case because the `ls` command failed each time (because I didn't have a file named `gefilte_fish`).

We can avoid the error message by using [`test`](http://pubs.opengroup.org/onlinepubs/9699919799/utilities/test.html#top).
Here `test -e` tests for the existence of a file.

```
$ test -e gefilte_fish && echo "we have fish"
$ test -e gefilte_fish || echo "we have no fish"
we have no fish
```

Here are some arguments to `test`:

* `-e pathname`: `pathname` is valid entry (directory, file, etc)
* `-n string`: string is non-empty
* `-z string`: string is empty
* `string1 = string2`: `string1` is identical to `string2`
* `string1 != string2`: `string1` is different than `string2`
* `i -eq j`: integer `i` is equal to `j`
* `i -ne j`: integer `i` is not equal to `j`
* `!` at the beginning negates what occurs to the right

The `!` works like so:

```
$ test ! -e gefilte_fish && echo "we have no fish"
we have no fish
```

You can use curly braces to delimit a compound statement:

```
test ! -e gefilte_fish && {
  echo "whitefish, pike, and carp"
  echo "are used to make gefilte fish"
}
```

You may be asking, *why would I do this, rather than use a regular if statement?*
Indeed, shell does have an if statement, which we'll cover next, and you *should* use that form when you are writing scripts, but I hope this helps you understand how these tools work together.


## [Control flow](http://pubs.opengroup.org/onlinepubs/9699919799/utilities/V3_chap02.html#tag_18_09_04)

We can rewrite our example above using an if statement:

```
if test -e gefilte_fish
then
  echo "we have fish"
else
  echo "we have no fish"
fi
```

I really want to emphasize that the "condition" of the if command (here a `test` statement) is a command, and we evaluate the subsequent statements based on the return status of that command.
After all, this is shell scripting, which is based on the results of running commands!

Now, to confuse things, shell introduces bracket notation, which is actually short-hand for a `test` call.
This makes it look like shell is a typical programming language:

```
if [ -e gefilte_fish ]
then
  echo "we have fish"
else
  echo "we have no fish"
fi
```

but in fact, `[` is just equivalent to calling `test` (if you don't believe me, try `man [`).
It's just more annoying, because you have to remember that trailing `]`.

There's also a while loop:

```
while test ! -e gefilte_fish
do
  echo "What, still no gefilte fish?"
  sleep 1
done
```

There's also a for loop.
Here's a silly example, which tells you about the file type of everything in your home directory.

```
for i in ~/*
do
  echo "what is the type of $i?"
  file $i
done
```

_Exercise:_ loop over the output of the `users` command on rhino, running `finger` on each user.


## The tiniest bit of `sed`

The use of `sed` is heresy.
If you need to do serious string manipulation, it's time for a general-purpose programming language.
But, if you don't mind hanging out in the pillory now and again, it can be handy in a pinch.

The command `sed` is short for stream editor.
It enables editing of streams of characters.
One related example that you might know about already is piping output of one command through `grep` to filter the results.

Here we look at processes happening on the system that contain the string `watchdog`.

```
ps aux | grep watchdog
```

Now, say we prefer cats over dogs.
We can replace every occurrence of dog with cat in this example like so:

```
ps aux | sed "s/dog/cat/" | grep watchcat
```

Protip: you can use any character for the delimiter of the regular expression, for example here `#`:

```
ps aux | sed "s#dog#cat#" | grep watchcat
```

which can be very handy if you want to replace strings that have slashes in them.



## Functions

Shell has functions.
The syntax is a little strange, in which the function is declared just by putting `()` after the name of the new function, and the variables are processed using special variables such as `$1`.

```
ldo () {
    $1 "$(ls -t | head -n 1)"
}
```

Actually, this one is pretty handy.
`$1` is the first argument.
So this function executes the first argument of the command on the result of running `ls -t | head -n 1`, which is the most recently touched file.
Might be handy for your `.bashrc` or `.zshrc` file.

As you might expect, `$2` is the second argument, `$3` the third, etc.
The variable `$@` is all of the arguments.

Just call a function like you would a regular shell command, e.g.

```
ldo head
```

will run `head` on the most recently modified file.

_Exercise:_ write a function, `lcd`, which enters the most recently touched directory. Note that you'll need to look at the command-line flags for `cd`.


## Command line arguments

Sometimes it's nice to be able to specify arguments to your script.
This is done just like for functions, with the same special variables `$1` and so on.
Note that `$#` is the number of arguments, scripts often start with something like:

```
if test $# -lt 2
then
    echo "usage: $0 fasta_file output_directory"
    exit
fi

fasta_file=$1
output_directory=$2
```

this checks that we have the right number of arguments and complains if we don't get them, then renames them something more reasonable.
Note that `$0` is the name of the script.

_Exercise:_ write a shell script that takes a directory as an argument, then loops at all of the files in that directory and unzips all of the zip files, gunzips all of the gzipped files, and leaves all other file types untouched.
