# Lecture 6: Working with data using tidyverse

Much of the time biologists spend while performing data analysis involves manipulating, reformatting, and transforming data. In this lecture, we'll use `tidyverse` functions from the `dplyr` package to perform these "data munging" tasks. We'll cover a few of the functions available for these tasks and learn a common programming structure (piping), which will equip you to continue developing powerful data manipulation code that will be foundational for visualization and statistical analysis.

## Learning objectives

- Select and apply functions from `tidyverse` to perform manipulations to tabular data
- Apply pipes as a programming structure to connect input and output from multiple functions

## Class materials

- Note that Github does not display `.html` files natively, so you need to download [`lecture.html`](lecture.html) and open in a file browser to view slides in presentation mode. [`lecture.org`](lecture.org) and [`reveal.js`](reveal.js) include the unformatted code and files used to create the slides; you should not need to use any of these files.
- RStudio has a [data manipulation cheatsheet](https://github.com/rstudio/cheatsheets/raw/master/data-transformation.pdf) that should help you identify what functions are useful for certain tasks.

## Reminders

- Homework 2 is currently available and is due Tuesday, October 22 at noon. You should have received an email containing an invitation to create your repository using GitHub Classroom. Contact Kate (khertwec at fredhutch.org) with any questions or concerns.
- We've had some questions about the difference between R Markdown documents and R Notebooks. The short answer is that R Notebooks are a specialized document written in R Markdown that updates as you run the code, while regular R Markdown documents need to be knit at the end of writing code for content to appear in the final document. More information can be found [here](http://uc-r.github.io/r_notebook). 
