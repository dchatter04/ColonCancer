# Gradient boosting reveals spatially diverse cholesterol gene signatures in colon cancer

**Authors:** Xiuxiu Yang<sup>1†</sup>, Debolina Chatterjee<sup>1†</sup>, Justin L Couetil<sup>2</sup>, Ziyu Liu<sup>3</sup>, Valerie D Ardon<sup>2</sup>, Chao Chen<sup>4</sup>, Jie Zhang<sup>2</sup>, Kun Huang<sup>1,2,5</sup>, Travis S Johnson<sup>1,5,6*</sup>


1. Department of Biostatistics and Health Data Science, Indiana University School of Medicine, Indianapolis, IN 46202, USA
2. Department of Medical and Molecular Genetics, Indiana University School of Medicine, Indianapolis, IN 46202, USA
3. Department of Statistics, Purdue University, West Lafayette, IN 47906, USA
4. Department of Biomedical Informatics, Stony Brook University, Stony Brook, NY11733, USA
5. Melvin and Bren Simon Comprehensive Cancer Center, Indianapolis University School of Medicine, Indianapolis, IN 46202, USA
6. Indiana Biosciences Research Institute, Indianapolis, IN 46202, USA


**Correspondence:** Travis S Johnson (johnstrs@iu.edu)

**Keywords:** colon cancer (CC), cholesterol, bile acids, prognostic genes, machine learning (ML), spatial transcriptomics (ST)

**Abstract:**

Colon cancer (CC) is the second most common cause of cancer deaths and the fourth most prevalent cancer in the United States. Recently cholesterol metabolism has been identified as a potential therapeutic avenue due to its consistent association with tumor treatment effects and overall prognosis. We conducted differential gene analysis and KEGG pathway analysis on paired tumor and adjacent-normal samples from the TCGA Colon Adenocarcinoma project, identifying that bile secretion was the only significantly downregulated pathway. To evaluate the relationship between cholesterol metabolism and CC prognosis, we used the genes from this pathway in several statistical models like Cox proportional Hazard (CPH), Random Forest (RF), Lasso Regression (LR), and the eXtreme Gradient Boosting (XGBoost) to identify the genes which contributed highly to the predictive ability of all models, ADCY5, and SLC2A1. We demonstrate that using cholesterol metabolism genes with XGBoost models improves the stratification of CC patients into low and high-risk groups compared with traditional CPH, RF, and LR models. Spatial transcriptomics (ST) revealed that SLC2A1 (glucose transporter 1, GLUT1) colocalized with small blood vessels. ADCY5 localized to stromal regions in both the ST and protein immunohistochemistry. Interestingly, both these significant genes are expressed in tissues other than the tumor itself, highlighting the complex interplay between the tumor and microenvironment, and that druggable targets may be found in the ability to modify how “normal” tissue interacts with tumors.

**Data and code availability:**
RNA-seq raw counts were retrieved from The Cancer Genome Atlas (TCGA) [TCGA data](https://portal.gdc.cancer.gov/projects) to study the relationship between cholesterol and CC prognosis. TCGA-COAD (N=512) [(Network, 2012)](https://doi.org/10.1038/nature11252)  read counts were normalized with the transcripts per million (TPM) method. After the data filtering process by removing the duplicates for each patient, 456 CC samples and 41 adjacent-normal tissues with survival information, age, gender, and stage were included for further analysis. The three Gene Expression Omnibus microarray datasets [GEO datasets](https://www.ncbi.nlm.nih.gov/geo/) were used for external validation cohorts: GSE17538 (N=232) [(Smith et al., 2010)](https://doi.org/10.1053/j.gastro.2009.11.005), GSE33113 (N=90) [(Felipe de Sousa et al., 2011)](https://doi.org/10.1016/j.stem.2011.10.008), and GSE39582 (N=566) [(Marisa et al., 2013)](10.1371/journal.pmed.1001453). The workflow for processing and analyzing the data is shown in Figure A below.

The preprocessed datasets used for the analysis in this paper can be found in the given Onedrive link [Datasets](https://indiana-my.sharepoint.com/:f:/g/personal/dchatter_iu_edu/Ev01OSB1ERlEpVfAN-MaJAEBN7wY4pUEQHnkkau6dTmqYg?e=tgZneM)

The R codes for the LASSO, Random Forest, and the XGBoost boost model can be found under the files list at [Rcodes](https://github.com/dchatter04/ColonCancer)

**Following is the workflow diagram for processing and analyzing the data.**

![Title](https://github.com/dchatter04/ColonCancer/blob/main/workflow.png)

