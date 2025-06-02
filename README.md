# Manifold Fitting Reveals Metabolomic Heterogeneity and Disease Associations in UK Biobank Populations
## Pre-requisites
* MATLAB R2022b
* R Version: 4.2.0 or higher
## Description
* ukb_nmr_main_process.R: The main process of our workflow.
* ukb_nmr_lifestyle_Cumulative_Incidence.R: Evaluations of the impact of sleep, physical activity, and smoking on the cumulative incidence of disease within a metabolically low-risk UK Biobank subgroup
* ukb_nmr_cox_diseases_build_Manifold.R: This function constructs an integrated dataset combining UMAP-based low-dimensional embeddings, clustering, and UK Biobank health and socioeconomic features to support downstream disease risk modeling such as Cox regression.
* ukb_nmr_disease_subgroups_cox.R: Cox regression across ICD-10 diseases.
## Citation
If you use this code for your research, please cite our papers.
```
@article{yao2024mfmetabolomicheterogeneity,
  title={Manifold Fitting Reveals Metabolomic Heterogeneity and Disease Associations in UK Biobank Populations},
  author={Li, Bingjie and Su, Jiaji and Lin, Runyu and Yau, Shing-Tung and Yao, Zhigang},
  journal={Manuscript},
  year={2024}
}
```
