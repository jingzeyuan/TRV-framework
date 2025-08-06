# TRV-framework

This repository accompanies the tutorial manuscript **"Detecting and Addressing Construct Bias in Cross-Cultural Psychological Measurement: A Tutorial"**, which introduces the Testing–Refining–Validating (TRV) framework. This framework provides a structured and theory-informed approach for identifying and addressing **structure-based construct bias** in cross-cultural psychological measurement.

Below is a description of the repository's contents and how each file supports the procedures introduced in the tutorial.

---

## 📁 Data Files

These three files contain the original survey data used in the tutorial's empirical demonstration:

- **China_Survey.xlsx**  
  Raw survey data collected from the Chinese sample.

- **India_Survey.csv**  
  Raw survey data collected from the Indian sample.

- **US_Survey.csv**  
  Raw survey data collected from the U.S. sample.

These datasets are merged and preprocessed into a single dataset (`Full_data`) within the main tutorial script.

---

## 🧠 Core Analysis Script

- **R codes for TVR framework_final.R**  
  This is the **main script** that implements the full TRV framework using the above datasets. It includes:
  
  - Data preparation and merging (creating `Full_data`)  
  - Testing phase: configural model  
  - Refining phase: structure modification using hACO  
  - Validating phase: model evaluation  
  - Real data analysis and annotated R code examples

> Refer to this script while working through the step-by-step tutorial in the manuscript.

---

## ⚙️ Custom Functions

Two additional R scripts contain reusable functions developed specifically for this framework. Their usage is explained and demonstrated in the tutorial:

- **hACO_MultiGroup_CFA.R**  
  This script includes the **hybrid Ant Colony Optimization (hACO)** algorithm designed to search for optimal model structures under the multigroup CFA context. It supports the Refining phase of the TRV framework.

- **Evaluate_invariance.R**  
  This script includes the `evaluate_invariance()` function used to assess measurement invariance across multiple groups during the Testing and Validating phases.

---

## 📄 README.md

This file. It provides an overview of the repository content and their relation to the tutorial.

---

## 🔗 Citation & License

If you use this code or framework in your research, please cite the accompanying tutorial paper (citation info will be updated after publication). All code is released for academic use under the [MIT License](https://opensource.org/licenses/MIT).

---

## 🧪 Reproducibility

All scripts are written in base R and depend on publicly available packages such as `lavaan`. The tutorial was tested under R 4.3.x. Please refer to the comments in each script for further guidance.

---

## 📬 Contact

For questions or collaborations, please contact [your email here].

