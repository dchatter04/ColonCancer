# ColonCancer

**Title of the paper:** Gradient boosting reveals spatially diverse cholesterol gene signatures in colon cancer


**Authors:** Xiuxiu Yang<sup>1†</sup>, Debolina Chatterjee<sup>1†</sup>, Justin L Couetil<sup>2</sup>, Ziyu Liu<sup>3</sup>, Valerie D Ardon<sup>2</sup>, Chao Chen<sup>4</sup>, Jie Zhang<sup>2</sup>, Kun Huang<sup>1,2,5</sup>, Travis S Johnson<sup>1,5,6*</sup>


1. Department of Biostatistics and Health Data Science, Indiana University School of Medicine, Indianapolis, IN 46202, USA
2. Department of Medical and Molecular Genetics, Indiana University School of Medicine, Indianapolis, IN 46202, USA
3. Department of Statistics, Purdue University, West Lafayette, IN 47906, USA
4. Department of Biomedical Informatics, Stony Brook University, Stony Brook, NY11733, USA
5. Melvin and Bren Simon Comprehensive Cancer Center, Indianapolis University School of Medicine, Indianapolis, IN 46202, USA
6. Indiana Biosciences Research Institute, Indianapolis, IN 46202, USA


**Correspondence:** Travis S Johnson (johnstrs@iu.edu)


**Keywords:** colon cancer (CC), cholesterol, bile acids, prognostic genes, machine learning (ML), spatial transcriptomics (ST)


**Abstract**

Colon cancer (CC) is the second most common cause of cancer deaths and the fourth most prevalent cancer in the United States. Recently cholesterol metabolism has been identified as a potential therapeutic avenue due to its consistent association with tumor treatment effects and overall prognosis. We conducted differential gene analysis and KEGG pathway analysis on paired tumor and adjacent-normal samples from the TCGA Colon Adenocarcinoma project, identifying that bile secretion was the only significantly downregulated pathway. To evaluate the relationship between cholesterol metabolism and CC prognosis, we used the genes from this pathway in several statistical models like Cox proportional Hazard (CPH), Random Forest (RF), Lasso Regression (LR), and the eXtreme Gradient Boosting (XGBoost) to identify the genes which contributed highly to the predictive ability of all models, ADCY5, and SLC2A1. We demonstrate that using cholesterol metabolism genes with XGBoost models improves stratification of CC patients into low and high-risk groups compared with traditional CPH, RF and LR models. Spatial transcriptomics (ST) revealed that SLC2A1 (glucose transporter 1, GLUT1) colocalized with small blood vessels. ADCY5 localized to stromal regions in both the ST and protein immunohistochemistry. Interestingly, both these significant genes are expressed in tissues other than the tumor itself, highlighting the complex interplay between the tumor and microenvironment, and that druggable targets may be found in the ability to modify how “normal” tissue interacts with tumors.


![Title](https://github.com/dchatter04/ColonCancer/blob/main/workflow.png "Title")

<img src="https://github.com/dchatter04/ColonCancer/blob/main/workflow.png" title="Title" width="500">
