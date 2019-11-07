# Lecture 13

For our next set of lectures, we'll be heading back to Python to explore some specific tools common in computational biology research. Lecture 13 and 14 from Jesse Bloom will introduce you to the use of packages in general, and [Biopython](https://biopython.org) in particular, which includes many useful tools for manipulation and analysis of data related to molecular biology.

## Learning objectives

After this section, you should be able to write Python functions that perform some useful real biological analyses:

- Get the reverse complement of a sequence (learn to use a Python `dict`).

- See if a sequence matches using ambiguous nucleotides (learn to use the Python `re` module).

- Parse barcodes from a real deep-sequencing experiment (learn to use `biopython`).

## Class materials

- The content for this lecture is contained in the Jupyter notebook [lecture13.ipynb](lecture13.ipynb) located in this directory.

- Please install [Biopython](https://anaconda.org/anaconda/biopython) prior to class by opening Terminal (Mac) or Anaconda Prompt (Windows) and executing the following code: `conda install -c anaconda biopython` . Alternatively (on either platform), open Anaconda Navigator, go to "Environments" and click "not installed". Search for biopython click the box to install.

- This class requires use of Python 3.6+. If you installed Anaconda recently (e.g., within the last few months, since the start of the course), you should be fine. If you are working with an older installation of Python, please see [these instructions](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-python.html#updating-or-upgrading-python) for checking your version and updating if necessary.

## Reminders

- Homeworks 4 and 5 were due at noon on Thursday, November 7.
