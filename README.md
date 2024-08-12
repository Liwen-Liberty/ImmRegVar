# Identify and quantify the effect of genomic alterations on CD8 T cell infiltration in tumor microenvironment based on interpretable deep learning models
===========================================================================


[![license](https://img.shields.io/badge/python_-3.11.5_-blue)](https://www.python.org/)
[![license](https://img.shields.io/badge/torch_-2.3.1_-blue)](https://pytorch.org/)
[![license](https://img.shields.io/badge/R_-4.2.2_-blue)](https://www.r-project.org/)

Tumor immunotherapy, particularly immune checkpoint inhibitors (ICI), has achieved revolutionary progress in cancer treatment, but challenges remain in predicting patient responses to therapy. To address this issue, we developed ImmRegVar, a Transformer-based model designed to predict CD8+ T cell infiltration by analyzing somatic genomic variations, including single nucleotide variants (SNVs) and copy number variations (CNVs), in 7972 solid tumor samples. Compared to traditional regression methods, ImmRegVar demonstrates significant superiority, achieving higher predictive accuracy and better model performance.More importantly, we introduced an innovative explanatory factor in our research. This factor skillfully integrates the attention weights from the Transformer model with feature ablation techniques, allowing for the precise quantification of the actual contribution of each gene mutation to the model’s predictions. This approach not only provides a clear and interpretable metric for understanding how specific gene mutations influence CD8+ T cell infiltration but also offers a new perspective on uncovering immune regulatory mechanisms. Additionally, through extensive random permutation tests, we validated the robustness and reliability of the explanatory factor, further confirming its stability in high-dimensional sparse data.Furthermore, the explanatory factor aligns with existing findings in medical literature regarding immune regulatory genes, enhancing its biological credibility. This makes our research not only theoretically innovative but also practically supportive for personalized immunotherapy. By focusing on the contribution of gene mutations, our study deepens the understanding of key factors in treatment, offering new pathways for researching immune regulatory mechanisms and their clinical applications. ImmRegVar and explanatory factor are freely available at https://github.com/Liwen-Liberty/ImmRegVar.

![Image text](https://github.com/Liwen-Liberty/ImmRegVar/blob/main/Figures/Figure1.png)

To identify and quantify the impact of genomic alterations on CD8 T-cell infiltration in the tumor microenvironment, we designed a neural network model based on the Transformer architecture. This model is intended to handle high-dimensional sparse genomic data and predict the effects of genomic alterations.First, we utilized a dataset containing gene identifiers and their corresponding mutation scores. Each gene is represented by a discrete identifier, while the mutation score indicates the extent of variation for each gene. Within this dataset, we integrated data on Single Nucleotide Variants (SNV) and Copy Number Variations (CNV).Then we built a neural network model based on the Transformer architecture.To enhance our understanding of genetic mutation data and improve the interpretability of model predictions, we have introduced a feature ablation method. In our research, each gene exists in two states: zero, indicating no mutation, and non-zero values representing various mutation scenarios. Due to the significantly higher frequency of non-mutation states compared to mutations, our data matrix forms a high-dimensional sparse matrix. In this context, the feature ablation method allows us to systematically assess the influence of each genetic feature on Transformer model predictions.

## Table of Contents

- [Installation](#installation)
- [Quick start](#quick-start)
- [Contributing](#contributing)
- [Cite](#cite)
- [Contacts](#contacts)
- [License](#license)


## Installation

ImmRegVar is tested to work under the:

```
* Python 3.11.5
* Torch 2.3.1
* R 4.2.2
* Numpy 1.24.3
* Other basic Python and r toolkits
```
### Installation of other dependencies
* Install [R package glmnetcr](https://github.com/cran/glmnetcr) using ` devtools::install_github("cran/glmnetcr") ` in the R environment if you encounter any issue.


# Quick start
To reproduce our results:


## 1. Application and verification
```
python ./main.py
```

**Arguments**:

| **Arguments** | **Detail** |
| --- | --- |
| **Solid_tumor_mutation_status_mat.csv** | Mutation status profile of 7972 solid tumor samples with matched expression data. |
| **Solid_tumor_CD8T_score_df.csv** | Cytolytic activity score of 7972 solid tumor samples. |


**Values**:

| **Output** | **Detail** |
| --- | --- |
| **model_state_dict.pt** | ImmRegVar model. |
| **attention_weights.npy** | attention weights matrix in ImmRegVar. |


## 2. Comparison of model explanatory factors and its robustness

See Python Code


Visualization of results:
<div align="center">
  <img src="https://github.com/Liwen-Liberty/ImmRegVar/blob/main/Figures/Figure2.png" alt="Editor" width="1000">
</div>



## 3. Important features selected based on explanatory factors

See Rscripts

Visualization of results:

<div align="center">
  <img src="https://github.com/Liwen-Liberty/ImmRegVar/blob/main/Figures/Figure3.png" alt="Editor" width="1000">
</div>


===========================================================================





# Contributing

yi He, liwen Xu, shaoliang Peng*

# Cite
<p align="center">
  <a href="https://clustrmaps.com/site/1bpq2">
     <img width="200"  src="https://clustrmaps.com/map_v2.png?cl=ffffff&w=268&t=m&d=4hIDPHzBcvyZcFn8iDMpEM-PyYTzzqGtngzRP7_HkNs" />
   </a>
</p>


# Contacts
If you have any questions or comments, please feel free to email: hyeliza4394@gmail.com.

# License

[MIT ? Richard McRichface.](../LICENSE)
