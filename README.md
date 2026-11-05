# nf1-single-cell-malignancy-classification
Reproducible scRNA-seq workflow for malignant cell identification and classification in NF1-associated tumors using machine learning.

# NF1 Single-Cell Malignancy Classification

**An R-to-Python computational workflow for single-cell transcriptomic analysis, biologically informed malignancy labeling, feature engineering, and machine learning in NF1-associated tumors.**

**Author:** Marta Morais Alonso  
**Master's Thesis:** Bioinformatics and Biostatistics — Universitat Oberta de Catalunya (UOC) / Universitat de Barcelona (UB)  
**Year:** 2026  
**Academic distinction:** Matrícula de Honor (highest academic distinction)

## Overview

Neurofibromatosis type 1 (NF1) is a genetic disorder associated with the development of peripheral nerve sheath tumors, including malignant peripheral nerve sheath tumors (MPNSTs). Distinguishing malignant from non-malignant cells is challenging because of tumor heterogeneity and the coexistence of multiple cellular populations.

This project develops a computational workflow to identify high-confidence malignant cells from single-cell RNA sequencing (scRNA-seq) data and evaluate supervised machine learning methods for their classification.

The work combines biological annotation, transcriptomic analysis, feature engineering, statistical evaluation, and robustness testing. Particular attention is given to label reliability and the risk of overestimating predictive performance when cells from the same biological sample are treated as independent observations.

The main technical contribution is the development of a structured downstream analysis workflow, from validation of the original `SingleCellExperiment` object to supervised model comparison and robustness assessment.

## Research objectives

- Validate and characterize a preprocessed single-cell transcriptomic dataset.
- Develop a conservative strategy for malignancy labeling using biological annotations and available evidence.
- Construct alternative transcriptomic feature representations for supervised learning.
- Compare Random Forest, linear SVM, and a neural network implemented in PyTorch.
- Assess model performance using both cell-level and sample-group-aware evaluation.
- Investigate model stability under controlled perturbations of the training data.

## Dataset

The project uses a previously processed and annotated `SingleCellExperiment` object supplied as the starting point of the academic research.

| Characteristic | Value |
|---|---|
| Cells | 117,133 |
| Genes | 9,131 |
| Tumor contexts | ANF, PNF, MPNST, and control nerve |
| Modality | scRNA-seq |
| Original format | RData / SingleCellExperiment |
| Python analysis format | AnnData / H5AD |

The input object already contained expression matrices, biological metadata, and dimensionality reductions. The project therefore focuses on downstream transcriptomic analysis rather than raw sequencing read processing.

## Computational methodology

### 1. R/Bioconductor exploration and data conversion

The original object was examined in R using `SingleCellExperiment` and `scater`.

This stage included:

- Inspection of expression assays, cell annotations, and dimensionality reductions.
- Exploratory t-SNE visualizations and cell-population summaries.
- Detection and correction of duplicated `colData` column names.
- Preparation of a working copy of the original object.
- A pilot conversion to AnnData/H5AD using `anndataR`.
- Validation of the exported data in Python with Scanpy.

The R scripts use project-relative paths with `here`. The Python analysis was developed in Google Colab using a project directory in Google Drive.

### 2. Conservative malignancy labeling

Three candidate labeling strategies were compared using tumor context, cellular annotations, `scna` information, and cluster-level evidence.

Instead of assigning a binary label to every cell, the final strategy distinguishes three categories:

| Cell category | Number of cells |
|---|---:|
| Malignant — high confidence | 20,555 |
| Non-malignant — high confidence | 19,889 |
| Uncertain | 76,689 |

This produces a supervised learning subset of **40,444 high-confidence cells**, while retaining uncertain cells outside the binary classification task.

The proposed labels were assessed through internal biological validation, including tumor and cell-type distributions, transcriptomic signatures, marker genes, and t-SNE representations.

### 3. Feature engineering

Three feature representations were constructed:

| Representation | Description | Features |
|---|---|---:|
| `all_genes` | Full log-transformed expression matrix | 9,131 |
| `scores` | Schwann, mesenchymal, and endothelial transcriptomic scores | 3 |
| `panel_plus_scores` | Selected marker-gene panel combined with transcriptomic scores | 19 |

The supervised model comparison focused on `scores` and `panel_plus_scores`.

Variables directly associated with the construction of malignancy labels were separated from the model predictors to reduce information leakage.

### 4. Machine learning

The project compares three classification approaches:

- **Random Forest**, providing a tree-based ensemble baseline.
- **Linear SVM**, implemented with standardized predictors.
- **Multilayer Perceptron (MLP)**, implemented and trained with PyTorch.

Two evaluation strategies were used:

- **Cell-level:** stratified train/test splits of individual cells.
- **Group-level:** repeated holdout evaluation using `GroupShuffleSplit`, with the `batch` identifier used to keep cells from the same group together.

The high-confidence dataset comprised 16 `batch`-defined groups. Five valid group-level splits were selected, requiring both classes to be represented in the training and test sets.

Group-level evaluation was prioritized because cells originating from the same biological sample may share biological and technical characteristics. This strategy reduces within-group information leakage, although the limited number of independent tumor samples constrains generalization.

Performance was evaluated using ROC-AUC, PR-AUC, accuracy, precision, recall, F1-score, and confusion matrices.

### 5. Robustness analysis

Additional experiments assessed the behavior of the three models under two types of perturbation:

- Reducing the size of the training set.
- Modifying the composition of non-malignant reference populations.

These experiments used the `panel_plus_scores` representation and group-level evaluation.

## Key results

The `panel_plus_scores` representation achieved better overall performance than the scores-only representation.

Under group-level evaluation, the MLP obtained the highest overall performance among the compared configurations.

| Model | ROC-AUC | PR-AUC | Accuracy | F1-score |
|---|---:|---:|---:|---:|
| Random Forest | 0.9919 | 0.9914 | 0.9360 | 0.9238 |
| Linear SVM | 0.9934 | 0.9929 | 0.9477 | 0.9407 |
| **MLP — PyTorch** | **0.9941** | **0.9936** | **0.9553** | **0.9511** |

*Mean results for `panel_plus_scores` across five valid group-level holdout splits.*

For the selected MLP configuration, the standard deviations across splits were **0.0061 for ROC-AUC** and **0.0505 for F1-score**. This variability, particularly in F1-score, should be considered when interpreting the results.

The findings support the feasibility of classifying high-confidence malignant and non-malignant cells within the analyzed dataset. They do not establish clinical diagnostic performance or generalization to independent cohorts.

## Project workflow

The main analytical sequence consists of 11 Python notebooks:

| Notebook | Main task |
|---|---|
| 01 | Pilot H5AD validation |
| 02 | Complete H5AD validation |
| 03 | Initial Scanpy exploration |
| 04 | Structural auditing and metadata dictionary |
| 05 | Malignancy labeling strategy |
| 06 | Biological validation of labels |
| 07 | Feature engineering and ML dataset construction |
| 08 | Random Forest and linear SVM baselines |
| 09 | PyTorch neural network |
| 10 | Comparative evaluation and reference model selection |
| 11 | Model stability under controlled perturbations |

The Python workflow is complemented by R scripts covering package installation, initial exploration, and preparation for H5AD export, as well as an R Markdown exploration report.

## Technology stack

**Languages:** Python, R

**Bioinformatics:** Bioconductor, SingleCellExperiment, scater, anndataR, AnnData, Scanpy

**Data analysis:** NumPy, pandas, SciPy, Matplotlib

**Machine learning:** scikit-learn, PyTorch

**Development environments:** RStudio, R Markdown, Jupyter Notebook, Google Colab

## Reproducibility and data availability

The research was completed using a combination of local RStudio workflows and Python notebooks executed through Google Colab.

The original transcriptomic dataset was supplied for the academic project, and its redistribution conditions have not yet been verified. The original RData file and derived H5AD datasets are therefore excluded from the public release.

The analytical workflow depends on these input files and on the original project directory structure. In the R export script, the pilot H5AD conversion is active, while the complete export is retained as an optional commented block.

The public repository is being prepared to document the code, analytical methodology, and results. It should not yet be considered a fully self-contained, independently executable pipeline.

## Limitations

The work is limited to the scRNA-seq modality; scATAC-seq integration was not implemented.

The malignancy labels are derived from conservative biological rules rather than independently established clinical ground truth. Consequently, the models are evaluated against operational labels constructed within the same study.

Although group-aware evaluation and perturbation experiments provide stronger internal evidence than cell-level evaluation alone, the limited number of independent tumor samples remains an important constraint, and independent external validation is necessary.

**This is a research methodology and computational classification study, not a clinically validated diagnostic tool.**

## Authorship

**Marta Morais Alonso — 2026**

Master's Thesis completed in June 2026, awarded Matrícula de Honor.

No software license has been assigned to this repository at this stage. The rights to third-party materials and the original biological dataset remain with their respective holders.
