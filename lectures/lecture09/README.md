# Lecture 9: First steps in the shell

**Erick Matsen ([@ematsen](https://twitter.com/ematsen), [matsen.group](http://matsen.group))**

Now that you have some experience using R for data manipulation and visualization, we'll move on to explore a few of the other common tools used in computational biology. This lesson introduces you to the command line (Unix shell), which is useful (and necessary) for automating tasks, working with files/directories, and using compute clusters. We'll be integrating our knowledge of the shell with Python coding starting with lecture 10, as these tools often work hand-in-hand.

## Learning objectives

After this course, you should be able to:

- Navigate a file system in the shell
- Run commands in the shell
- Understand basic usage of piping and redirection
- To be able to edit a file using vim
- To write a basic shell script
- To be able to run a command in parallel using GNU Parallel

## Class materials

- All instructions for this lesson are available in this repository.
  - See the Tutorial section below to get started with the exercises
  - View the introductory presentation by opening the [`slides.html`](slides/slides.html) in a web browser from your local copy of the repository 
- This material requires use of a [unix shell](https://github.com/fredhutchio/tfcb_2019/tree/master/software#unix-command-line). Different "flavors" of shell have slight variations in commands available. These materials are designed to work on a compute cluster at Fred Hutch called rhino. Please see [these instructions](https://github.com/fredhutchio/tfcb_2019/tree/master/software/unix_rhino.md) for logging on to rhino, and note there is an extra step to log in off campus. You can execute most of these commands on your own computer (e.g., without logging in to rhino), but don't be surprised if some of the commands and options don't work!

## Reminders

- Homework 3 (genomic data in R) is available through GitHub Classroom and is due Tuesday, October 28 at noon.

## Tutorial

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


This is an interactive tutorial to teach you about navigating the shell.
For reference you might be interested in the excellent [course material put out by Software Carpentry](https://swcarpentry.github.io/shell-novice/).

To get started, make sure you're connected to the Marconi wireless network. Open Unix shell program on your computer and execute the following commands (where username is your HutchNetID):

    ssh username@rhino
    git clone https://github.com/fredhutchio/tfcb_2019.git
    cd tfcb_2019/lectures/lecture09

Now you can start the first lesson by clicking on the `first-steps.md` link in the file list above and starting your adventure!

If you have problems at any point, flag us down and we'll come by to help out.
