# RNA Expression Analysis Tool

A Python-based web application for exploring and analyzing RNA expression datasets. The application allows users to upload RNA expression data, identify biological groups within the dataset, and visualize expression patterns interactively.

## Features
Upload RNA expression datasets in .tsv, .csv, or .txt format
Process and organize gene-expression data using Pandas
Identify and compare biological sample groups
Visualize RNA expression patterns across samples
Interactive web interface built with Streamlit
Designed as a foundation for future biomarker and differential-expression analysis

## Technologies
Python
Pandas — data processing and analysis
Streamlit — interactive web application
Machine Learning — pattern analysis within biological datasets
Matplotlib / visualization tools — expression visualization

## Dataset

The repository includes an RNA expression dataset from the NCBI Gene Expression Omnibus (GEO), allowing the application to be tested using real biological data.

## Project Goal

The goal of this project is to develop a computational tool for exploring gene-expression data and identifying potentially meaningful patterns between biological groups. Future development will focus on differential-expression analysis, biomarker identification, and pathway-level analysis.

## Running Locally

Clone the repository and install the required dependencies:

pip install -r requirements.txt

Run the application with:

streamlit run app.py

The application will open in a local web browser.

Future Development
Differential-expression analysis
Gene-specific expression lookup
Statistical significance testing
Biomarker identification
Pathway enrichment analysis
Support for additional RNA-seq datasets
