# 1. deTS introduction (v1.0)
Genome-wide association studies (GWAS) and next-generation sequencing technologies have identified hundreds of thousands of disease-associated variants and genes. Interpretation of these variants could be greatly enhanced in tissue-specific systems. However, there are many diseases or traits where the causal tissues or cell types remain unknown. deTS: Tissue-Specific Enrichment Analysis is an R package to identify the most relevant tissues for candidate genes or for gene expression profiles. deTS builds on two pre-processed reference panels. We implemented different statistic tests for different forms of query data. 
# 2. Usage
## 2.1 Installing deTS
### Requirements
deTS relies on R (>= 3.4), pheatmap (>= 1.0.10), RColorBrewer (>= 1.1)  
The pheatmap relies on CRAN. Please follow their installation instruction.  
`> install.packages("pheatmap")  `
### To download the codes, please do:
`git clone https://github.com/bsml320/deTS.git`  
`cd deTS`  
Then open the R:   
`R`  
`> install.packages("deTS_1.0.tar.gz")  `
### deTS loading
Load the deTS package and dependent library  
`> library(deTS)`  
`> library(pheatmap)`  
## 2.2 Built-in data loading
deTS requires two reference panels to conduct the enrichment test: one from GTEx and the other from ENCODE. For GTEx, a matrix including the summary statistics for each tissue is also needed. All datasets have been included in the package. After installation of the package, one can load the data using the following commands:  

Load the t-statistic matrix for the GTEx panel:  
`> data(GTEx_t_score)`  
Load the z-score matrix for the ENCODE panel:  
`> data(ENCODE_z_score)`  
Then `"GTEx_t_score"` and `"ENCODE_z_score"` will be loaded to R enviroment.  
## 2.3 Input data
deTS deals with two types of enrichment analysis for different forms of query data. For convenience, we provide two Tissue-Specific Enrichment Analysis (TSEA) functions for query gene lists (single sample and multiple samples), and another function for RNA-Seq expression profiles tissue-specific enrichment analysis.    
### 2.3.1 TSEA for gene lists
When the query data are lists of genes, the Fisher’s Exact Test is implemented. The function is `tsea.analysis()`. The input is a vector of gene symbols. Here we used disease-associated genes identified from GWAS summary statistics as an example. The gene symbols can be found here:  
Load gene symbol from deTS package:  
`> data(GWAS_gene)`  
`> query.genes = GWAS_gene`  
Or you can read your own gene symbol list from a text file:  
`> dat = read.table("data/Gene_list.txt", head = F)`  
`> query.genes_user = dat[,1]`  

Nextly, we perform tissue-specific enrichment analysis for query gene list:  
`> tsea_t = tsea.analysis(query.genes, GTEx_t_score, ratio = 0.05, p.adjust.method = "bonferroni")`  
Here, the ratio is a value to define tissue-specific genes (default is 5%) and provides the first way of categorizing genes.  
The second way of grouping genes is based on the query genes. The two ways of category form a two by two table, which is used in the Fisher’s Exact Text.  
The Fisher's Exact Test results between query gene list and each tissue specific genes will be stored in variable `tsea_t`.  
You can check tissue-specific enrichment analysis result by:    
`> head(tsea_t)`  
For better visualization and summary, we provide one plot and one summary function to list the top 3 enriched tissues, simply run:  
`> tsea.plot(tsea_t, threshold = 0.05)`  
`> tsea_t_summary = tsea.summary(tsea_t)`  

### 2.3.2 TSEA for multiple gene lists  
In most condition, you might want to analysis multiple samples together, then you can upload a 0~1 table. In the table, gene labeled with 1 indicated significant associate within a sample, while 0 indicated not in a given sample. You can check the format of example data.  
Load multiple gene symbol from deTS package:  
`> data(GWAS_gene_multiple)`  
`> query.gene.list = GWAS_gene_multiple`  
Or you can read your own gene symbol list from a text file:  
`> dat = read.table("data/Gene_list_multiple.txt", head = T, row.names = 1)`  
`> query.gene.list_user = dat`  

To keep result reliable, please keep at least 20 genes for each samples.   
You can check the total genes number for each sample:  
`> colSums(query.gene.list)`  

Then, we can make tissue specific enrichment analysis for multiple samples by `tsea.analysis.multiple()` and plot the result by `tsea.plot()`. You can summary the top 3 most associated tissues by `tsea.summary()` function and save your result in to a text-format spreadsheet:  
Tissue-specific enrichment analysis in GTEx panel:  
`> tsea_t_multi = tsea.analysis.multiple(query.gene.list, GTEx_t_score, ratio = 0.05, p.adjust.method = "BH")`  
Save tissue-specific enrichment analysis result:  
`> write.csv(tsea_t_multi,"GWAS_multi_TSEA_in_GTEx_panel.csv")`  
Save the tissue-specific enrichment analysis plot:  
`> pdf ("GWAS_multi_TSEA_in_GTEx_panel.pdf", 6, 6, onefile = FALSE)`  
`> tsea.plot(tsea_t_multi, threshold = 0.05)`  
`> dev.off()`   
Save summary result in to a spreadsheet:  
`> tsea_t_multi_summary = tsea.summary(tsea_t_multi)`  
`> write.csv(tsea_t_multi_summary,"GWAS_multi_summary_GTEx_panel.csv")`

### 2.3.3 TSEA for RNA-Seq profiles
For a quick start, user can use ENCODE example RNA-Seq profiles make TSEA in GTEx panel:  
Load ENCODE query data:  
`> data(query_ENCODE)`  
`> query.matrix = query_ENCODE`  
Load correction variable:  
`> data(correction_factor)` 

As RNA-Seq samples are often heterogeneous, before in-depth analysis, it’s necessary to decode tissue heterogeneity to avoid samples with confounding effects. However, the raw discrete RPKM value should be normalized to continuous variable meet the normal distribution before t-test. We provided two normalization approaches: `"z-score"` and `"abundance"` in function `tsea.expression.normalization()`:  
(1) `"z-score"` normalization will calculate a z-score for the query sample for each tissue in the reference panel as below: e_i=(e_0-μ_t))/sd_t, where μ_t and sd_t were the mean and SD of tissue t.   
(2) `"abundance"` normalization will provide an abundance correction approach for the query sample for each tissue in the reference panel as below: e_i=(log2(e_0+1)/(log2(u_t+1)+1).  

We have the preloaded the test RPKM variable in `query.matrix` and correction variable in `correction_factor`, we take `"abundance"` normalization approach as an example, simply type:  
RNA-Seq profiles scale by abundance normalization:  
`> query_mat_abundance_nor = tsea.expression.normalization(query.matrix, correction_factor, normalization = "abundance")`  

After get normalized RPKM value, we submit it for `tsea.expression.decode()`:   
`> tseaed_in_GTEx = tsea.expression.decode(query_mat_abundance_nor, GTEx_t_score, ratio = 0.05, p.adjust.method = "BH")`  
Then, the tissue specific enrichment analysis for query RNA-Seq is finish. After tissue specific enrichment decode analysis, one-side t-test results between query RNA-Seq sample tissue specific genes (top 5%) versus remains genes (95%) is stored in variable `tseaed_in_GTEx`. Further analysis for top 3 most associated tissues is similar to previous analysis:  
`> tsea.plot(tseaed_in_GTEx, threshold = 0.05)`  
`> tseaed_in_GTEx_summary = tsea.summary(tseaed_in_GTEx)`  
`> write.csv(tseaed_in_GTEx_summary,"RNAseq_summary_in_GTEx_panel.csv")`  

To prove the robustness of our proposed pipeline, user can validate the two reference panels through self-validation. Simply, load GTEx example RNA-Seq profiles and perform tissue-specific enrichment analysis in ENCODE panel:  
`> data(query_GTEx)`  
`> query.matrix = query_GTEx`  
RNA expression profiles z-score normalization:   
`> query_mat_zscore_nor = tsea.expression.normalization(query.matrix, correction_factor, normalization = "z-score")`  
RNA expression profiles TSEA in ENCODE panel:  
`> tseaed_in_ENCODE = tsea.expression.decode(query_mat_zscore_nor, ENCODE_z_score, ratio = 0.05, p.adjust.method = "BH")`  

The reader is encouraged to open and view the file in a spreadsheet software, or inspect it directly within R using the command `fix(tseaed_in_ENCODE)`. In addition, sometime, you might want to edit some parameters for your own data, e.g., you can change the `GTEx_t_score` to `ENCODE_z_score` for ENCODE tissue specific enrichment analysis, you can also change the tissue specific genes `ratio` from `0.05` to `0.2`, or change the `p.adjust.method` to `"bonferroni"`.  

Further analysis for top 3 most associated tissues is similar to previous analysis:  
`> tsea.plot(tseaed_in_ENCODE, threshold = 0.05)`  
`> tseaed_in_ENCODE_summary = tsea.summary(tseaed_in_ENCODE)`  
`> write.csv(tseaed_in_ENCODE_summary,"RNAseq_summary_in_ENCODE_panel.csv")`  

## Citation
https://www.ncbi.nlm.nih.gov/pubmed/30824912

Pei G., Dai Y., Zhao Z., Jia P. (2019) deTS: Tissue-Specific Enrichment Analysis to decode tissue specificity. Bioinformatics, 35(19):3842-3845.
