# Homework 8: Capstone

## Objectives

This assignment is designed to help you:
1. Plan and implement a project in computational biology from start to finish
2. Assess and use previously published biological research data
3. Apply reproducible computational tools to summarize and visualize data to answer research questions
4. Document project progress using version control and a README file

## Summary

For your capstone project, you will select a dataset and identify a series of questions to answer using these data. You will apply one or more of the coding approaches covered this semester (or others you choose to learn on your own) to answer these questions. Your final submission will be a GitHub repository (`homework08`, from GitHub Classroom) containing code, figures, and documentation. 

Your performance on this assignment will be assessed using the attached [rubric](#rubric) (100 points total). Your assignment is due on Tuesday, December 10 at noon.

## Details 

To complete this project, you should:
1. Choose one of the following datasets. These datasets are all previously published and available through Data Dryad, a public repository. Please see the next section on identifying research questions for additional context in selecting a dataset.
- Glaucoma diagnosis https://datadryad.org/stash/dataset/doi:10.5061/dryad.q6ft5 
- Triple negative breast cancer imaging https://datadryad.org/stash/dataset/doi:10.5061/dryad.32765 
- Immunohistochemical typing of adenocarcinomas https://datadryad.org/stash/dataset/doi:10.5061/dryad.g8h71 
2. Identify three questions to answer using tools applied this semester. Appropriate questions include performing data filtering/assessment/grouping, summary statistics and assessing patterns using data visualization. Please note:
- The data you select from the list above may not follow best practices in data organization. You may need to convert file formats of the data prior to analyzing. We recommend you assess (and document) things like column headers and encoding of missing data as a part of your analysis.
- Some of these datasets are not particularly well-documented. The original manuscripts may provide additional context, but you will have the opportunity to describe the assumptions you made about the data and what could've been done to improve their reusability.
- These data come from previously published research studies. You are not required to ask the same questions as the original publications; this assignment is deliberately open-ended to allow you freedom to select from the tools we've used this semester.
- Your questions should require more than a minimal level of coding effort to assess (e.g., it is unacceptable to ask “What is the mean of [some column] in the dataset?”). However, you are also not required to apply advanced statistical testing in your analysis beyond what we have explicitly covered in class (though this is a useful skill). 
- You should include at least two data visualizations in your project, including application of best practices in visualizing as described in the first part of our course.
- An (oversimplified) dataset with questions/results: If you chose a dataset representing student enrollment in 200 different courses throughout ten years, your questions would be 1) Are all classes present across the entirety of the timeframe? (result is a yes or no, if no, include the number of classes that aren’t present), 2) What is the average enrollment in the largest 10 classes across ten years? (result is a table and line graph), 3) Does enrollment in each course vary with total enrollment across all classes? (result is a density plot with overlay of total enrollment)
3. Set up a project repository in GitHub. Your repository should:
- Be created from accepting the link from GitHub Classroom for homework08. 
- Represent appropriate practices for a project directory (subdirectories, file naming, etc) as described in the first few classes.
- Not include large data files (>50 MB). You will not be able to push large data files to GitHub. We recommend use of a `.gitignore` file, storing large files in other locations outside your project directory, and use of smaller/filtered/reduced intermediate data files when appropriate. 
- Track a reasonable version history of your project as it develops. In other words, you should not have only an initial commit and then upload of the final files.
4. Document your data, code, and other project components with a `README` file. This should follow best practices covered in the first part of the semester, including:
- Explanations of the directory structure and all included files.
- Description of the data (with citation) and where it can be found online
5. Write code to answer your three research questions. Your code should include:
- A markdown (`report.md` or `report.Rmd`) or ipython notebook (`report.ipynb`) serving as a report that includes the following:
  - A meaningful title and brief (markdown) introduction to the report.
  - A written (markdown) section at the top titled “About the data” that describes any assumptions you made about the data and what would’ve improved your ability to reuse it.
  - Each research question clearly stated, associated code, results/conclusions/interpretations (which may include print statements, markdown text chunks, and data visualizations).
  - Complete code comments that describes the tasks being accomplished.
  - For each visualization, a justification for the type of visualization and your formatting choices.
  - A written (markdown) section at the end titled “Reproducibility” that comments on the ease of reproducibility of your analysis and the analysis in the original paper, which reflects on the concepts covered in the first few weeks of class.
-  In addition to your report, you may (but are not required to) include additional files used for analysis (either as intermediate results, or scripts). Make sure they are documented in your README!
        
## Rubric
