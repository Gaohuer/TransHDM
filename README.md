# TransHDM

**High-Dimensional Mediation Analysis via Transfer Learning**

`TransHDM` is an R package that provides a framework for **high-dimensional mediation analysis using transfer learning**. The package integrates large-scale **source datasets** to enhance the detection power of potential mediators in **small-sample target studies**.

The proposed method addresses **data heterogeneity** across studies through transfer regularization and debiased estimation, while maintaining control of the false discovery rate.

The methodology implemented in TransHDM is based on the framework proposed in Pan et al. (2025) <doi:10.1093/bib/bbaf460>.

---

# Installation

You can install the development version of `TransHDM` from GitHub using the **devtools** package.

```r
install.packages("devtools")
devtools::install_github("Gaohuer/TransHDM")
```

Alternatively, using **remotes**:

```r
install.packages("remotes")
remotes::install_github("Gaohuer/TransHDM")
```

---

# Functions

### Main Function

* `TransHDM()`

The main function of the package performs **transfer learning–based high-dimensional mediation analysis** by integrating information from source datasets to improve mediator detection in the target dataset.

### Data Generation

* `gen_simData_homo()`
  Generate **homogeneous simulation datasets** for mediation analysis.

* `gen_simData_hetero()`
  Generate **heterogeneous simulation datasets** with distributional differences between source and target data.

### Diagnostics

* `source_detection()`
  Detects **reliable source datasets** and evaluates **source-target heterogeneity**.

### Baseline Methods

* `lasso()`
  Implements the **LASSO-based mediator selection** method.

* `dblasso()`
  Implements the **debiased LASSO approach** for mediation analysis.

### Variable Screening

* `SIS()`
  Performs **Sure Independence Screening (SIS)** for dimension reduction in high-dimensional settings.

# Examples

Detailed usage examples, including simulation studies and analysis workflows, are provided in the **package vignettes**.

After installing the package, you can access them using:

```r
browseVignettes("TransHDM")
```

---

# Citation

If you use **TransHDM** in your research, please cite the associated article:


> Pan L, Liu Y, Huang C, Lin R, Yu Y, Qin G. Transfer learning reveals the mediating mechanisms of cross-ethnic lipid metabolic pathways in the association between APOE gene and Alzheimer's disease. Brief Bioinform. 2025;26(5):bbaf460. <doi:10.1093/bib/bbaf460>


You can also obtain the citation information directly from R:

```r
citation("TransHDM")
```

---

# License

This package is released under the **GPL-3 License**.
