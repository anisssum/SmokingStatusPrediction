# A predictive model for predicting the smoking status

**Author**  
- Anna Chesnokova

**Supervisor**
- Ivan Valiev

## Table of contents

- [Introduction](#introduction)
- [Pipeline](#pipeline)
  - [Data preparation for differential gene expression analysis](#data-preparation-for-differential-gene-expression-analysis)
  - [Differential gene expression analysis](#differential-gene-expression-analysis)
  - [Building a logistic regression model](#building-a-logistic-regression-model)
  - [Integration of mutation signature data](#integration-of-mutation-signature-data)
  - [Building machine learning models](#building-machine-learning-models)
- [Results](#results)
- [Summary](#summary)
- [Reference](#reference)

## Introduction

Tobacco Smoking History Status:
- Lifelong Non-smoker (less than 100 cigarettes smoked in Lifetime) = 1
- Current smoker (includes daily smokers and non-daily smokers or occasional smokers) = 2
- Current reformed smoker for > 15 years (greater than 15 years) = 3
- Current reformed smoker for ≤15 years (less than or equal to 15 years) = 4
- Current reformed smoker, duration not specified = 5

Lung cancer remains one of the most prevalent and deadly forms of cancer worldwide [^1].
Smoking is a major risk factor, influencing both the initiation and progression of lung cancer,
creating certain patterns of mutations and gene expression [^2].
Recent advancements in genomics and transcriptomics have provided valuable insights into the
molecular mechanisms underlying smoking-related carcinogenesis [^3][^4].
This study aims to develop a machine learning model to predict the smoking status
of lung cancer patients based on their genetic and transcriptomic profiles.
The primary goal is to enhance the accuracy of classification by integrating transcription factors (TFs),
differential gene expression data (DEGs) and mutational signatures.

## Pipeline

Gene expression, smoking status data and mutation annotation files were retrieved from the database UCSC XENA [^5].

In order to use `.ipnyb` notebooks you will need to download ['requirements.txt'](https://github.com/anisssum/SmokingIndexPrediction/blob/main/requirements.txt).

### Data preparation for differential gene expression analysis

1. Obtaining columns with tumor type and smoking status from ['tcga.tcga.LUAD.MD'](https://github.com/anisssum/SmokingIndexPrediction/blob/main/data/tcga.tcga.LUAD.MD):
['transpose_ID_SmokingHistory.csv'](https://github.com/anisssum/SmokingIndexPrediction/blob/main/data/transpose_ID_SampleType.csv) and
['transpose_ID_SmokingHistory.csv'](https://github.com/anisssum/SmokingIndexPrediction/blob/main/data/transpose_ID_SmokingHistory.csv).

2. For the read count [tcga.gene_sums.LUAD.R109.gz](https://github.com/anisssum/SmokingIndexPrediction/blob/main/data/tcga.gene_sums.LUAD.R109.gz) (forced compressed for github),
a file was obtained to translate the gene IDs [gene_annotation_table.txt](https://github.com/anisssum/SmokingIndexPrediction/blob/main/data/gene_annotation_table.txt)
from the gene annotation file [human.gene_sums.R109.gtf](https://github.com/anisssum/SmokingIndexPrediction/blob/main/data/human.gene_sums.R109.gtf).

### Differential gene expression analysis

Script for data preparation and calculation of differential gene expression: [DESeq_LUAD.R](https://github.com/anisssum/SmokingIndexPrediction/blob/main/DESeq_LUAD.R).

1. Select only primery and recurent tumors.

2. Prepare Metadata for DESeq2 [^6] and create a table for machine learning models ([data_t_sort_all.txt.gz](https://github.com/anisssum/SmokingIndexPrediction/blob/main/data/data_t_sort_all.txt.gz), (forced compressed for github)).

3. Obtaining train id ([train_id.csv](https://github.com/anisssum/SmokingIndexPrediction/blob/main/data/train_id.csv)) for calculating differential gene expression.

4. Create the DESeq2 dataset using the counts matrix and metadata.

5. Normalize the data using size factors.

6. Saving DESeq2 results: [DESeq_res.csv](https://github.com/anisssum/SmokingIndexPrediction/blob/main/data/DESeq_res.csv).

7. Save the top of 250 differentially expressed genes to CSV files (padj < 0.05):
[DEG_up.csv](https://github.com/anisssum/SmokingIndexPrediction/blob/main/data/DEG_down.csv), [DEG_up.csv](https://github.com/anisssum/SmokingIndexPrediction/blob/main/data/DEG_up.csv).

### Building a logistic regression model

Construction of a logistic regression model using the expression of transcription factors and DEGs as an example, similar was done with other chipsets: [DF_TF_logreg.ipynb](https://github.com/anisssum/SmokingIndexPrediction/blob/main/DF_TF_logreg.ipynb).

1. Prepare data and obtain train IDs for DESeq2.

2. Count tpm for genes with the realized function `read_counts2tpm()`,
that uses gene lengths ([gene_length_1.txt](https://github.com/anisssum/SmokingIndexPrediction/blob/main/data/gene_length_1.txt)) obtained from the gene annotation file
([human.gene_sums.R109.gtf](https://github.com/anisssum/SmokingIndexPrediction/blob/main/data/human.gene_sums.R109.gtf)) using a script [gene_length.bash](https://github.com/anisssum/SmokingIndexPrediction/blob/main/gene_length.bash).

3. Save the table with tpm value: [X_tpm_DF.csv](https://github.com/anisssum/SmokingIndexPrediction/blob/main/data/X_tpm_DF.csv).

4. Build logistic regression models for DEGs (tpm value and log10+1 tpm value).

6. Build logistic regression models for TFs obtained from the Enrichr database (tpm value and log10+1 tpm value) [^7].

### Integration of mutation signature data

The section of work with mutation annotation files is located in the folder `/maf`.

1. In the original file ([TCGA-LUAD.mutect2_snv.tsv](https://github.com/anisssum/SmokingIndexPrediction/blob/main/maf/data/TCGA-LUAD.mutect2_snv.tsv))  [^8]
select only single-nucleotide variants ([TCGA-LUAD.mutect2_sbs.tsv](https://github.com/anisssum/SmokingIndexPrediction/blob/main/maf/data/TCGA-LUAD.mutect2_sbs.tsv))
using a script [filter_sbs.sh](https://github.com/anisssum/SmokingIndexPrediction/blob/main/maf/filter_sbs.sh).

2. With a script [count_signatures.R](https://github.com/anisssum/SmokingIndexPrediction/blob/main/maf/count_signatures.R)
using the library `deconstructSigs` [^9] we will create a three nucleotide matrix and calculate the weight of each signature
[combined_table_signatures_tobacco.tsv](https://github.com/anisssum/SmokingIndexPrediction/blob/main/maf/data/combined_table_signatures_tobacco.tsv).

### Building machine learning models

Сonstruction of models using the expression of DEGs with the weight of each signature as an example, similar was done with other chipsets: [DEG_SIGS_models.ipynb](https://github.com/anisssum/SmokingIndexPrediction/blob/main/maf/DEG_SIGS_models.ipynb).

1. Replace patient IDs and combine data from DEGs and mutation-signature weights.

2. Build machine learning: Logistic Regression, Decision Tree Classifier, Random Forest Classifier, and LGBM Classifier [^10].

3. Calculating shap values.

## Results

### PCAplots

![PCAplot1](https://github.com/anisssum/SmokingIndexPrediction/blob/main/img/PCAplot_2.png)

_Figure 1. PCAplot on normalized read counts._

![PCAplot2](https://github.com/anisssum/SmokingIndexPrediction/blob/main/img/PCAplot_1.png)

_Figure 2. PCAplot by normalized log2 of the number of reads._

### Volcano plot for DEGs

![Volcano_Plot](https://github.com/anisssum/SmokingStatusPrediction/blob/main/img/plot_zoom_png.png)

### Transcription factors (ARCHS4):

- RHOXF2
  
- GATA1
  
- CTCFL
  
- MAEL

- TFDP3
  
- ZNF595

### F1 value for different feature sets

![F1_value](https://github.com/anisssum/SmokingStatusPrediction/blob/main/img/2024-05-27_01-35-50.png)

### SHAP value (effect on model result)

![SHAP_value](https://github.com/anisssum/SmokingStatusPrediction/blob/main/img/2024-05-27_01-49-05.png)

## Summary

The results indicate that DEGs and mutational signatures are powerful predictors of smoking history in lung cancer patients.
TFs, while biologically significant, did not provide sufficient discriminatory power on their own, likely due to the complexity of transcriptional regulation and the indirect nature of TF activity.
The superior performance of the combined feature set underscores the importance of integrating multiple data types to capture the multifaceted impact of smoking on lung cancer biology.
Future studies should focus on validating these findings in larger independent cohorts, improving machine learning models and exploring the potential for clinical application in personalized medicine.

## Reference:

[^1]:	Siegel, R. L., Miller, K. D., and Jemal, A. "Cancer statistics, 2020." CA: A Cancer Journal for Clinicians, 2020. doi:10.3322/caac.21590.

[^2]: Govindan R. et al. Genomic landscape of non-small cell lung cancer in smokers and never-smokers //Cell. – 2012. – Т. 150. – №. 6. – С. 1121-1134, doi:10.1016/j.cell.2012.08.024.

[^3]: Alexandrov, L. B., Nik-Zainal, S., Wedge, D. C., et al. "Signatures of mutational processes in human cancer." Nature, 2013. doi:10.1038/nature12477.

[^4]: Goldman, M.J., Craft, B., Hastie, M. et al. Visualizing and interpreting cancer genomics data via the Xena platform. Nat Biotechnol (2020), https://doi.org/10.1038/s41587-020-0546-8.

[^5]: Love MI, Huber W, Anders S (2014). “Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.” Genome Biology, 15, 550, doi:10.1186/s13059-014-0550-8.

[^6]: Tate J. G. et al. COSMIC: the catalogue of somatic mutations in cancer //Nucleic acids research. – 2019. – Т. 47. – №. D1. – С. D941-D947, https://doi.org/10.1093/nar/gky1015.

[^7]: Kuleshov, M. V., Jones, M. R., Rouillard, A. D., et al. "Enrichr: a comprehensive gene set enrichment analysis web server 2016 update." Nucleic Acids Research, 2016. doi:10.1093/nar/gkw377.

[^8]: Rosenthal, R., McGranahan, N., Herrero, J., et al. "DeconstructSigs: delineating mutational processes in single tumors distinguishes DNA repair deficiencies and patterns of carcinoma evolution." Genome Biology, 2016. doi:10.1186/s13059-016-0893-4.

[^9]: Pedregosa, F., Varoquaux, G., Gramfort, A., et al. "Scikit-learn: Machine learning in Python." Journal of Machine Learning Research, 2011.

[^10]: Lambert S. A. et al. The human transcription factors //Cell. – 2018. – Т. 172. – №. 4. – С. 650-665, https://doi.org/10.1016/j.cell.2018.01.029.
