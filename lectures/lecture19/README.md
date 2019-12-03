# Lecture 19: Putting it all together

This last session will walk through a realistic example of  analyzing high-throughput sequencing data that will use several of the tools that we learned in the previous 18 lectures.

## Learning objectives

After this session, you should be able to:
- outline a workflow for analyzing data from a new experiment.
- write `Python` scripts to count barcodes from high-througput sequencing reads.
- write `R` scripts to plot the summary count statistics.
- track your code in a `GitHub` repository.

## Class materials

- [`lecture.html`](lecture.html) contains the slides that we will start with.
- Note that Github does not display `.html` files natively, so you need to download [`lecture.html`](lecture.html) and open in a file browser to view slides in presentation mode. [`lecture.org`](lecture.org) and [`reveal.js`](reveal.js) include the unformatted code and files used to create the slides; you should not need to use any of these files.
- You should log into Rhino using the instructions from [lecture 9](https://github.com/fredhutchio/tfcb_2019/tree/master/lectures/lecture09#tutorial) for this session.
- The example workflow described here is adapted from the analysis carried out in [Park 2019](https://github.com/rasilab/ribosome_collisions_yeast).
