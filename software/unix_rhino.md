# Logging on to rhino

## On Fred Hutch's campus

Check to make sure you're connected to the Marconi wireless network. Open the program you'll be using to access Unix on your computer (Terminal on Mac or WSL/Git bash on Windows) and type the following (where `username` is your HutchNetID):

    ssh username@rhino

Hit enter to execute the command. The first time you execute this on a computer, you'll receive a response similar to:

> The authenticity of host 'rhino.fhcrc.org' can't be established.
> Are you sure you want to continue connecting (yes/no)?

Type `y` and press enter.

You'll be prompted to enter your username and then password (note that your cursor will not move as you type your password; this is normal!).

After you have successfully logged in, you'll see the information about the cluster printed to your screen. You'll be ready to start entering commands when you see a prompt like the following appear:

> username@rhino2:~$


## Off campus log-in

For more information about remote login, please see [this entry](https://sciwiki.fredhutch.org/scicomputing/access_methods/#access-via-a-remote-location) in the Fred Hutch Biomedical Data Science Wiki.

The short version: logging in off campus requires an additional step to connect to the campus network (where username is your HutchNetID):

    ssh username@snail.fhcrc.org

You'll see a message printed to the screen that starts with:

>  Welcome to the Fred Hutchinson Cancer Research Center

Then you can enter:

    ssh rhino

You'll then be able to interact with the cluster as if you were on campus.
