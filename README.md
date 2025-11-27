# DCDM Group 2: IMPC Phenotype Database & Explorer

## Overview

This repository contains the complete data management pipeline for the **DCDM Group 2** coursework. The project integrates high-throughput knockout mouse phenotyping data from the [International Mouse Phenotyping Consortium (IMPC)](https://www.mousephenotype.org/) with human disease ontologies.

The solution consists of three main components:

1.  **ETL Pipeline:** A series of R scripts to clean, normalise, and categorise raw phenotypic data.
2.  **Relational Database:** A 3NF-normalised MySQL schema connecting Genes, Procedures, Parameters, Analyses, and Diseases.
3.  **Interactive Dashboard:** An R Shiny application for exploring genotype-phenotype associations, performing dimensionality reduction (PCA), and visualising gene clusters.

## Repository Structure

```text
DCDM_Group2/
├── config/                 # Configuration files for file paths
│   └── paths.yaml
├── database2/              # Database artifacts and SQL dumps
│   └── DUMP                # MySQL database dump file
├── outputs/                # Generated CSVs for SQL import and intermediate data
├── scripts/                # Source code for Data Cleaning, SQL Prep, and Shiny App
│   ├── app.R                                   # Main R Shiny Dashboard application
│   ├── IMPC_APP.R                              # (Backup/Alternative) App script
│   ├── R_script_for _SQL_import.Rmd            # Final preparation of CSVs for SQL LOAD DATA
│   ├── STEP1_rawdata_establish.Rmd             # Initial raw data loading
│   ├── STEP2_cleaning_Disease_information.Rmd  # Cleaning disease ontology data
│   ├── STEP2_cleaning_IMPC_Procedure.Rmd       # Cleaning procedure metadata
│   ├── STEP2_cleaning_parameter_description.Rmd# Text mining parameter descriptions
│   └── STEP3_parameter_groupings.Rmd           # Categorising parameters (e.g., "Metabolism")
├── SMOKE_TEST.txt          # System integrity test file
└── .gitignore              # Git ignore configuration

## Prerequisites

* **R & RStudio:** (Version 4.0+)
* **MySQL Server:** (Version 8.0+)
* **R Packages:**
    * `tidyverse` (Data manipulation)
    * `shiny` (Web framework)
    * `DT` (Interactive tables)
    * `plotly` (Interactive charts)
    * `viridis` (Colour scales)
    * `readr`, `yaml` (File IO)

## Setup Instructions

### 1. Data Processing (ETL)
The scripts in the `scripts/` directory should be run in the following order to reproduce the dataset:

1.  **Raw Data Setup:** Run `STEP1_rawdata_establish.Rmd` to load the source files.
2.  **Cleaning Modules:** Run the three `STEP2` scripts to clean Procedures, Disease info, and Parameter descriptions.
3.  **Categorisation:** Run `STEP3_parameter_groupings.Rmd` to generate the logical groupings for the phenotypes.
4.  **SQL Preparation:** Run `R_script_for _SQL_import.Rmd`. This generates the final normalised CSV files (e.g., `Genes_for_SQL.csv`, `Analyses_for_SQL.csv`) in the `outputs/` folder.

### 2. Database Initialisation
1.  Create the database `DCDM_IMPC_Group2` in MySQL.
2.  Execute the table creation SQL commands (found inside `R_script_for _SQL_import.Rmd` or the SQL dump).
3.  Import the CSV files generated in the previous step using `LOAD DATA INFILE` commands.
    * *Note: Ensure `local_infile` is enabled if loading from client.*

### 3. Launching the Dashboard
To explore the data:
1.  Open `scripts/app.R` in RStudio.
2.  Ensure the input CSVs (`IMPC_cleaned3.csv`, `IMPC_parameter_groupings.csv`) are accessible to the script.
3.  Click **"Run App"**.

## Features

* **Gene View:** Univariate jitter plots showing significance scores ($-\log_{10} P$) across phenotype categories.
* **Phenotype View:** Comparative horizontal bar charts ranking genes by effect size for specific parameters.
* **Clustering:** Hierarchical clustering heatmaps identifying pleiotropic gene signatures.
* **PCA:** Principal Component Analysis to visualise gene separation based on phenotypic profiles.

## Authors
**DCDM Group 2**
