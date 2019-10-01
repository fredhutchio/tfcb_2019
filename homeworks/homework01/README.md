# Homework 1

This homework will assess your ability to use GitHub, organize projects, and apply tidy data principles, and represents content covered in lectures 01-03. You will also have the opportunity to familiarize yourself with GitHub Classroom, through which you will submit homework assignments for the duration of the course.

At the end of this assignment, you will have two public repositories. One repository (`tfcb-homework01`) will be a revision of `messy-project-directory/`, following the instructions below, and will be located in your own GitHub account. The second repository (prefixed with `homework01`) will belong to the GitHub organization set up specifically for homework collection; this repository will be automatically set up after you accept the invitation (from clicking the link sent via email). The second repository will contain a separate copy of this file that you will modify for your homework submission. Your answers will replace the lines specified in _italics_ below. Your answers in that file will be recorded automatically at the time due (October 8 at noon) for your homework submission.

## Problem 1

**10 points**

Make a GitHub account and populate your bio by including your research interests, place of work, location, professional webpage (or science-focused twitter account!). Here's an example [github.com/trvrb/](https://github.com/trvrb/).

_Include the link to your GitHub profile here._

## Problem 2

**10 points**

This question assesses your ability to track a project using Git.

Download the course files with the following link: [github.com/fredhutchio/tfcb_2019/archive/master.zip](https://github.com/fredhutchio/tfcb_2019/archive/master.zip). Copy the contents of `messy-project-directory/` to a new directory on your computer named `tfcb-homework01`.

Open GitHub Desktop and make a "New repository" with name `tfcb-homework01`. Set "Local Path" to the location in your computer where your `tfcb-homework01` directory can be found.

Create an initial commit that adds all the local files in `tfcb-homework01`:
- `Survey Data.xlsx`
- `get Species_list.py`
- etc...

You do not need to submit anything for this problem; your success will be evaluated in the next problem.

## Problem 3

**10 points**

This question assesses your ability to publish projects to GitHub.

Publish `tfcb-homework01` to GitHub using GitHub Desktop. Make sure to set this to be a "public" repository. The resulting repository can now be accessed at github.com/{your_name}/tfcb-homework01. If necessary, you can make this repository public by going to "Settings" from this page.

_Include the link to your `tfcb-homework01` GitHub repository here._

## Problem 4

**10 points**

This question assesses your ability to organize files and directories associated with research projects.

Organize files into a more consistent structure. Group images into a `images/` directory. Separate source code and data. Rename files to remove spaces and improve consistency.

Commit changes and publish to your public GitHub repository. Locate the URL of this commit by clicking on "commits" from your project page on GitHub, which should be similar in format to: https://github.com/fredhutchio/tfcb_2019/commit/16b5235bd2d908c96f22297813e6aaf9f172ad41

_Include the commit URL for your reorganized project here._

## Problem 4

**10 points**

This question assesses your ability to write a README with markdown formatting.

Create a file called `README.md` and populate with Markdown. Demonstrate headers, lists, links, embedded images (by linking to images contained in the directory) and tables in this readme.

Commit this file and publish to your public GitHub repository.

_Include the link to your `README.md` here._

## Problem 5

**10 points**

This question assesses your understanding of tidy data principles.

Clean up the file that was originally named `Survey Data.xlsx`. Some points to remember: aim for a single tidy data frame in a single tab, don't use formatting as data, use preferred date format, properly record null values.

Export this as a tab-delimited `.tsv` text file with Unix line endings.

Commit the modified `.xlsx` file and the `.tsv` file and publish to your public GitHub repository.

_Include the link to your `.tsv` file on GitHub here._
