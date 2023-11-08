# Single-cell-transcriptomic-molecular-classification-of-breast-cancer

#Background
Three clinical subtypes of breast cancer: Luminal (ER+, PR+/−), HER2+(HER2+, ER+/−, PR+/−) and triple negative (TNBC; ER −, PR −, HER2 −)
Breast cancer PAM50 molecular subtypes: batch transcriptomic analysis of PAM50 gene markers classifies breast cancer into five "intrinsic" molecular subtypes: lumen like (LumA and LumB), HER2 rich, Basal like and Normal like
There is~70-80% consistency between molecular subtypes and clinical subtypes. Although PAM50 provides important insights for prognosis and treatment, due to the limitations of batch transcriptomics, it is not possible to determine the subtype of breast cancer specific tumor cells, so it is currently necessary to develop a single cell breast cancer molecular typing method.

![PAM50](result/PAM50.png)

#Obtain breast cancer tumor cells
After quality control of single cell samples, published gene markers were used to further confirm cell annotation, and all major cell types were representative in all tumors and clinical subtypes. Visualize epithelial cells clearly separated from tumors through UMAP. Since breast cancer is mainly driven by DNA copy number changes, InferCNV was used to estimate the single cell copy number variation (CNV) spectrum to distinguish tumor cells from normal epithelial cells. Select T cells and endothelial cells as reference genes (Figure 2, 4, 8, 10, 18), as these types of cells rarely exhibit copy number variation.

![cell_annotation](result/cell_annotation.png)

![infercnv](result/infercnv.png)

