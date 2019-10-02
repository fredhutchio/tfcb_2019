# Lecture 3: Introduction to Data

**Trevor Bedford ([@trvrb](https://twitter.com/trvrb), [bedford.io](https://bedford.io))**

## Learning objectives

- Apply computational project organization principles specifically to good practices in data management
- Implement best practices for organizing tabular (spreadsheet-style) data according to tidy data principles

## Class materials

1. [File organization and naming](#project-and-data-organization)
2. [Tidy data](#tidy-data)

This class requires Microsoft Excel (or an equivalent program that can open `.xlsx` files); see [Software installation](../../software/) for more information.

## Reminders

- Recommended reading is available throughout each section
- Please have R and RStudio installed before the next class session. See [Software installation](../../software/) for more information. **Please also follow the instructions to [install tidyverse](../../software/README.md#installing-tidyverse).**
- Your first homework assignment is currently available and is due Tuesday, October 8 at noon. You should have received an email containing an invitation to create your repository using GitHub Classroom. Contact Kate (khertwec at fredhutch.org) with any questions or concerns.

## Project and data organization

It's important to keep a tidy project directory, even if something is not as the stage of being versioned on GitHub.

Some general advice:

1. **Encapsulate everything within one directory, which is named after the project.** Have a single directory for a project, containing all of the data, code, and results for that project. This makes it easier to find things, or to zip it all up and hand it off to someone else.
2. **Separate the data from the code.** I prefer to put code and data in separate subdirectories. I’ll often have a `data/` subdirectory and a `scripts/` (or `src`) subdirectory.
3. **Use relative paths (never absolute paths).** If you encapsulate all data and code within a single project directory, then you can refer to data files with relative paths (e.g., `../data/some_file.csv`). If you were to use an absolute path (like `~/Projects/SomeProject/data/some_file.csv` or `C:\Users\SomeOne\Projects\SomeProject\data\some_file.csv)` then anyone who wanted to reproduce your results but had the project placed in some other location would have to go in and edit all of those directory/file names.
4. **Write dates as YYYY-MM-DD**. This sorts properly and also avoids ambiguities.
5. **Include readme files**. This bit of documentation greatly helps in describing what a folder contains.
6. **Continually up-to-date directory.** I aim for a clean and up-to-date directory that is continually modified rather than a chronological directory structure (as described in [this article](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000424)). A separate electronic lab notebook with chronological entries can be hugely helpful for record-keeping purposes.

### File names

Borrowing excellent slide deck from Ciera Martinez and colleagues: [Reproducible Science Workshop: File Naming](https://rawgit.com/Reproducible-Science-Curriculum/rr-organization1/master/organization-01-slides.html#1)

### File organization

Continuing slide deck from Ciera Martinez and colleagues: [Reproducible Science Workshop: Organization](https://rawgit.com/Reproducible-Science-Curriculum/rr-organization1/master/organization-02-slides.html)

### Documenting your data

Ideally, your `data/` directory will include an additional `README` that (at bare minimum) includes a description of the data (e.g., what the rows and columns represent). Fully documented metadata (data about the data) will include:
- **Data description**: descriptions for labels of rows and columns (e.g., abbreviations defined, units of measurement)
- **Intent**: project description and purpose, reference to other published works and associated documentation (e.g., manuscripts)
- **Collection method**: experimental design
- **File structure**: if multiple data files are present, the organization and relationship among these files
- **Data type**: file types and formats
- **Provenance**: how data were created and by whom (including analytical pipelines), what normalization or transformations are currently represented

Documenting data can be a time-consuming process, but is often required to submit data to repositories. Since data publishing is a requirement for most academic research as a part of publication, keeping track of this information early on can save you time later, and increase the chances of other researchers using your data later (which means more citations for you).

### Miscellaneous advice

[More excellent advice from Karl Broman](https://kbroman.org/dataorg/)

## Tidy data

Tidy data is term from Hadley Wickham and refers to:

>A standard method of displaying a multivariate set of data is in the form of a data matrix in which rows correspond to sample individuals and columns to variables, so that the entry in the ith row and jth column gives the value of the jth variate as measured or observed on the ith individual.

Data in this form is much easier to deal with programmatically. This is also known as a _data frame_. This [tutorial presents a nice overview](https://r4ds.had.co.nz/tidy-data.html).

Observations as rows and variables as columns is an excellent standard to adhere to.

1. Each variable forms a column
2. Each observation forms a row.
3. Each type of observational unit forms a table

See for example, single cell RNA sequencing data, with cells as rows and genes as columns. This is also the way that relational databases (MySQL, Postgres, etc...) are constructed.

### Exercise on tidy data

1. Demonstrate conversion of simple example dataset. Work from [Table 2 in Bedford et al. 2014](https://bedford.io/papers/bedford-flux/), available as an [Excel table in the course repo](tables/influenza-evolutionary-parameters.xlsx).

2. Pair up to work from an HI (haemagglutination-inhibition) table and convert to tidy data. Data available as an [Excel table in the course repo](tables/influenza-titer-data.xlsx).

## File formats

Saving data as plain text files is necessarily to process this data with either R or Python. You can export from Excel to `.tsv` (tab-delimited, my preferred format) or `.csv` (comma-delimited). A few things to note when exporting data files in these formats:
- Beware that [line endings](https://support.nesi.org.nz/hc/en-gb/articles/218032857-Converting-from-Windows-style-to-UNIX-style-line-endings) differ between Windows and Unix (including Mac), though the text editors we recommend for this class can deal with this
- Exporting from Excel only works for the currently displayed spreadsheet. If you have multiple sheets, you'll need to export multiple times.

## More detailed examples

* [Return to Sidney's dengue example](https://github.com/blab/dengue-antigenic-dynamics)

### Further reading

Some suggested readings include:

* [Good Enough Practices for Scientific Computing by Wilson et al.](https://swcarpentry.github.io/good-enough-practices-in-scientific-computing/) as already referenced
