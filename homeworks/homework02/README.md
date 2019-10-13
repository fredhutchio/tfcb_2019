TFCB 2019: Homework 2
=====================
Oct 15, 2018

*Due 12pm, Oct 22, 2018*

In this homework, we will learn how to analyze a recently published deep sequencing dataset using `tidyverse` functions.

In the process, we will learn some new functions in `tidyverse` and apply them to our data analysis.

You can open the [homework02.Rmd](homework02.Rmd) file in this directory to write your answers inside `RStudio` and push the knitted
`.md` and associated image files to Github Classroom.

-   [Problem 1](#problem-1)
-   [Problem 2](#problem-2)
-   [Problem 3](#problem-3)
-   [Problem 4](#problem-4)
-   [Problem 5](#problem-5)
-   [Problem 6](#problem-6)
-   [Problem 7](#problem-7)
-   [Problem 8](#problem-8)
-   [Problem 9](#problem-9)
-   [Problem 10](#problem-10)


``` r
library(tidyverse)
```

Problem 1
---------

**10 points**

Provide a &lt;100 char description and a URL reference for each of these functions.

1.  `!`
2.  `is.na`
3.  `is.numeric`
4.  `anti_join`
5.  `desc`
6.  `slice`
7.  `all_vars`
8.  `funs`
9.  `filter_if`
10. `mutate_if`

Problem 2
---------

**10 points**

Add a comment above each code line below explaining what the code line does or why that code line is necessary.

Keep each comment to less than 2 lines per line of code and &lt; 80 chars per line.

``` r
annotations <- read_tsv("ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/locus_groups/protein-coding_gene.txt") %>% 
  select(ensembl_gene_id, symbol, name, gene_family, ccds_id) %>% 
  filter(!is.na(ccds_id)) %>% 
  print()
```

``` r
data <- read_tsv("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE89nnn/GSE89183/suppl/GSE89183_Counts.txt.gz") %>% 
  rename(ensembl_gene_id = `ENSEMBL gene`) %>%
  print()
```


Problem 3
---------

**10 points**

1.  Convert the plot below to use `log10` instead of linear scales for both axes.
2.  Show axis tick labels as 10<sup>0</sup>, 10<sup>1</sup>, 10<sup>2</sup>,10<sup>3</sup>, 10<sup>4</sup>, 10<sup>5</sup> for each axis.
3.  There are two many points overlapping in certain regions. Use a different `geom_` function to convey to your reader how many overlapping points are present in each region.

``` r
data %>% 
  select(CD34_shRPL5_RNA_1, CD34_shRPS19_RNA_1) %>% 
  ggplot(aes(x = CD34_shRPL5_RNA_1, y = CD34_shRPS19_RNA_1)) +
  geom_point()
```

![](README_files/figure-markdown_github/unnamed-chunk-4-1.png)

Problem 4
---------

In problems 4--6, after each operation, assign the result back to the `data` variable.

**10 points**

Select the following columns from the `data` variable you created above. `ensembl_gene_id`, `CD34_shRPL5_RPF_1`, `CD34_shRPL5_RPF_2`, `CD34_shRPL5_RNA_1`, `CD34_shRPL5_RNA_2`, `CD34_shRPS19_RPF_1`, `CD34_shRPS19_RPF_2`, `CD34_shRPS19_RNA_1`, `CD34_shRPS19_RNA_2`, `CD34_shLuc_RPF_1`, `CD34_shLuc_RPF_2`, `CD34_shLuc_RNA_1`, `CD34_shLuc_RNA_2`.

Problem 5
---------

**10 points**

Filter the result from Problem 4 to include only rows where each of the 12 numerical columns you selected above has 50 counts or more. This is a simple way to avoid genes that have very low counts. You might be tempted to do this step separately for each of the 12 columns. But you can write lot less code if you use the `filter_if` and `all_vars` functions you learned above.

Problem 6
---------

**10 points**

After filtering in Problem 5, divide each of the above 12 numerical columns by the corresponding median value in each column. This median normalization is typically done in high-throughput experiments to normalize for sample-to-sample difference in read depth. Again, you can write lot less code if you use the `mutate_if` and `funs` functions you learned above.

Problem 7
---------

**10 points**

After we do the above filtering and median-normalization, let us calculate translation efficiency as the average ratio of the RPF and RNA reads for each treatment condition. Then we calculate how this translation efficiency changes between target (`rpl5` and `rps19`) and control (`luc`) shRNAs.

The code implementing the above steps is shown below, but it has a few errors. Correct them.

``` r
lfc <- data %>% 
  mutate(mean_rpl5_te = ((CD34_shRPL5_RPF_1 + CD34_shRPL5_RPF_2) / 
                            (CD34_shRPL5_RNA_1 + CD34_shRPL5_RNA_2)) %>% 
  mutate(mean_rps19_te = ((CD34_shRPS19_RPF_1 + CD34_shRPS19_RPF_2) / 
                            (CD34_shRPS19_RNA_1 + CD34_shRPS19_RNA_2)) %>% 
  mutate(mean_shluc_te = ((CD34_shLuc_RPF_1 + CD34_shLuc_RPF_2) / 
                            (CD34_shLuc_RNA_1 + CD34_shLuc_RNA_2)) %>% 
  select(ensembl_gene_id, mean_rpl5_te, mean_rps19_te) %>% 
  mutate(lfc_te_rpl5 == log2(mean_rpl5_te / mean_shluc_te),
         lfc_te_rps19 == log2(mean_rps19_te / mean_shluc_te)) %>% 
  print()
```

Problem 8
---------

**10 points**

Create a new dataframe called `mean_lfc` from `lfc` containing a new column called `avg_lfc`. `avg_lfc` should be the average of the log2 fold-change in TE (`lfc_te`) upon knockdown of RPL5 and RPS19.

Then select only the gene id column and the new column that you just created.

Problem 9
---------

**10 points**

Join the above `mean_lfc` dataframe with the `annotations` dataframe created at the top of the document.

Problem 10
----------

**5 points**

Select only the bottom 10 genes with the lowest `avg_lfc` and display the gene `symbol`, gene `name` and `avg_lfc` for these genes.

5 bonus points if you use one of the new functions you learned above!
