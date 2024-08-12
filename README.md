# Identify and quantify the effect of genomic alterations on CD8 T cell infiltration in tumor microenvironment based on interpretable deep learning models
===========================================================================


[![license](https://img.shields.io/badge/python_-3.11.5_-blue)](https://www.python.org/)
[![license](https://img.shields.io/badge/torch_-2.3.1_-blue)](https://pytorch.org/)
[![license](https://img.shields.io/badge/R_-4.2.2_-blue)](https://www.r-project.org/)

Tumor immunotherapy, particularly immune checkpoint inhibitors (ICI), has achieved revolutionary progress in cancer treatment, but challenges remain in predicting patient responses to therapy. To address this issue, we developed ImmRegVar, a Transformer-based model designed to predict CD8+ T cell infiltration by analyzing somatic genomic variations, including single nucleotide variants (SNVs) and copy number variations (CNVs), in 7972 solid tumor samples. Compared to traditional regression methods, ImmRegVar demonstrates significant superiority, achieving higher predictive accuracy and better model performance.More importantly, we introduced an innovative explanatory factor in our research. This factor skillfully integrates the attention weights from the Transformer model with feature ablation techniques, allowing for the precise quantification of the actual contribution of each gene mutation to the model’s predictions. This approach not only provides a clear and interpretable metric for understanding how specific gene mutations influence CD8+ T cell infiltration but also offers a new perspective on uncovering immune regulatory mechanisms. Additionally, through extensive random permutation tests, we validated the robustness and reliability of the explanatory factor, further confirming its stability in high-dimensional sparse data.Furthermore, the explanatory factor aligns with existing findings in medical literature regarding immune regulatory genes, enhancing its biological credibility. This makes our research not only theoretically innovative but also practically supportive for personalized immunotherapy. By focusing on the contribution of gene mutations, our study deepens the understanding of key factors in treatment, offering new pathways for researching immune regulatory mechanisms and their clinical applications. ImmRegVar and explanatory factor are freely available at https://github.com/Liwen-Liberty/ImmRegVar.

![Image text](https://github.com/Liwen-Liberty/ImmRegVar/blob/main/Figures/Figure1.png)
