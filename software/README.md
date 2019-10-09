# Software installation

We ask you come prepared to class with a laptop on which you can participate in coding activities (talk to Kate if you don't have access to a computer and we can arrange one to borrow). Please follow the instructions below to install the required software for this course (all programs should come pre-installed on your computer, or are freely available for academic use). We'll use the other files in this directory to test the software together in class. The tools you'll need include:

- [Text Editor](#text-editor) for simple file viewing / manipulation
- [Spreadsheet program](#spreadsheet-program), like Microsoft Excel
- [Git and GitHub](#git) for version control (please note you'll also need to share your GitHub username with us!)
- [R and RStudio](#r) for R statistical programming
- [WSL or Terminal](#unix-command-line) for Unix command line
- [Jupyter Notebooks via Anaconda](#python) for Python

## Text Editor

Biological data is almost exclusively represented as text, and we will be writing code in text files. It's useful to be able to open, examine, and edit text files using a light-weight text editor. Microsoft Word is not suitable for these types of files, and we do not recommend the text editors set as default on your computer, as these programs are often not optimized for working with code. If you have a text editor already installed that you have previously used for code, you are welcome to use that software. If not, please install one of our recommendations below, which are appropriate for viewing very large text files and doing find-and-replace tasks using regular expressions.

### MacOSX

BBEdit can be downloaded from [here](https://www.barebones.com/products/bbedit/download.html) or installed via the App Store.

### Windows

[Notepad++](https://notepad-plus-plus.org) can be downloaded by going to the [Download page](https://notepad-plus-plus.org/downloads/) and selecting the latest release.

## Spreadsheet program

Spreadsheet programs are a useful way for us as humans to interact with data. The most common of these is Microsoft Excel. Commands may differ a bit between programs, but the general ideas for thinking about spreadsheets are the same. If you are working on a computer owned by Fred Hutch, Microsoft Office (including Excel) is available through the Self Service application. If you are working on a personal computer that doesn’t have a spreadsheet program, you can use a free, open source program called LibreOffice.

Install LibreOffice by going to the installation page. The version for your operating system should automatically be selected. Click Download Version X.X.X (whichever is the most recent version). You will go to a page that asks about a donation, but you don’t need to make one. Your download should begin automatically.
Once the installer is downloaded, double click on it and LibreOffice should install.

## Git

Git is version control software, which helps you keep track of changes made to files. GitHub is a repository for data and code tracked with Git, and is a mechanism for publishing and collaborating on project development.

### GitHub Account

If you do not already have one, please register for a [GitHub account](https://github.com). Please note that your name and email will be publicly visible through GitHub by default, but more information on controlling privacy settings can be found [here](https://help.github.com/en/articles/setting-your-commit-email-address).

**Please share your GitHub username with us [here](https://docs.google.com/forms/d/e/1FAIpQLSdtlewc-Py9BZ33_Kgi-u2d2ilHYEDEYniwKVjW8mIp508r8A/viewform?usp=sf_link) so we can add you to GitHub classroom (used for homework submission).**

### GitHub Desktop App

The website for [GitHub Desktop](https://desktop.github.com) should auto-detect your operating system and allow you to download and install the software. If you are using a Windows machine, please select the option to install command line tools.

## R

R and RStudio are separate downloads.
R is the "engine", while RStudio is an integrated desktop environment (IDE) that makes using R much more pleasant.
R must be installed before RStudio.
Follow the instructions below for your operating system to install them.
If you are working on a computer owned by Fred Hutch,
RStudio + R is available through the Self Service application.

### Windows

- Download the installer for the latest version of R from [CRAN](http://cran.r-project.org/bin/windows/base/release.htm).
  The file will begin downloading automatically.
- Double-click the downloaded `.exe` file and follow the prompts to install.
- Go to the [RStudio download page](https://www.rstudio.com/products/rstudio/download/#download).
- Under _Installers_, click the link for the _Windows Vista/7/8/10_ installer to download it.
- Double-click the downloaded `.exe` file and follow the prompts to install (default options are acceptable).
- Once both are installed, launch RStudio and make sure there are no error messages.

### MacOSX

- Download the installer for the latest version of R compatible with your version of macOS from [CRAN](https://cran.r-project.org/bin/macosx/).
  If you are not using a recent version of macOS you may have to scroll down to _Binaries for legacy OS X systems_ and find the one appropriate for your version of macOS.
  To check what version of macOS you are using, click the apple icon in the upper left corner of your screen and go to _About This Mac_.
  Please note the instructions on that page for downloading and installing [XQuartz](https://www.xquartz.org/) if necessary.
- Double-click the downloaded `.pkg` file and follow the prompts to install (default options are acceptable).
- Go the the [RStudio download page](https://www.rstudio.com/products/rstudio/download/#download).
- Under _Installers_, click the link for the your OSX version's installer to download it.
- Double-click the downloaded `.dmg` file, then open the RStudio folder that appears on your desktop. Drag the RStudio icon into the Applications folder.
- Once everything is installed, launch RStudio and make sure there are no error messages.

### Installing tidyverse

- Open Rstudio. 
- Click the Packages tab in the lower right panel. 
- Click the Install button (upper left corner of the panel). In the empty space for Packages, type `tidyverse`. The other defaults (Install from CRAN and the Install to Library path) should be ok. Make sure the box next to "Install dependencies" is checked, and click Install.
- If your installation is successful, you should see tidyverse appear in the list below. 

## Unix command line

### Windows

Windows 10 comes with a new feature called Windows Subsystem for Linux (WSL) that allows you to access Unix tools on your computer. Installation instructions can be found [here](https://docs.microsoft.com/en-us/windows/wsl/install-win10).

Another option (such as if you are not running Windows 10) is Git for Windows, which also installs Git command-line tools. You can download [here](https://git-scm.com/download/win) and install with default options.

### MacOSX

Macintosh operating systems are built on Unix, so many of the tools you'll need are pre-installed on your computer. You can access the command line through an application called **Terminal**. You can either search for this in Finder, or use the Go drop-down menu to locate it in the Utilities folder.

## Python

We will use [Jupyter notebooks](http://jupyter.org) to record code, output, and text throughout the course.
We recommend installing Python using Anaconda,
which includes Jupyter notebooks and most of the other packages we'll use for the course, according to the following instructions:
- Download the [Anaconda](https://www.anaconda.com/download/) installer for
Python 3.x for your particular operating system.
- Double-click the downloaded file and follow the prompts to install Anaconda (default options are acceptable).
- For assistance troubleshooting installation, please go [here](https://jupyter.readthedocs.io/en/latest/install.html).
