# Whole-Blood Transcriptomics Data-Based Disease Classification Pipeline

## Project Description
This repository contains the code and results for a pipeline designed for disease classification based on whole-blood transcriptomics data. The pipeline aims to process transcriptomics data, perform analysis, and generate classification results and visualizations.

## Contents
-   `codes/`: Contains R scripts for the pipeline functions and the main execution script.
    -   `pipeline_function.R`: Core functions used in the data processing and classification pipeline.
    -   `WBT_DC_pipeline.R`: The main script to run the whole-blood transcriptomics disease classification pipeline.
-   `results/`: Stores the output of the pipeline, including ROC comparison plots and other analysis results.
    -   `roc_comparison_*.pdf`: PDF files containing Receiver Operating Characteristic (ROC) curve comparisons for different disease classifications.

## Setup and Installation
To run this pipeline, you will need R and several R packages.
 **Clone the repository:**
    ```bash
    git clone https://github.com/Mengzhensw/WBT-DC-Pipeline-Whole-Blood-Transcriptomics-data-based-Disease-Classification.git
    cd WBT-DC-Pipeline-Whole-Blood-Transcriptomics-data-based-Disease-Classification
    ```

## Usage
1.  **Prepare your data:** All processed data required to run the WBT-DC pipeline have been deposited in Zenodo (https://zenodo.org/records/17990588)
2.  **Run the pipeline:**
    Execute the main pipeline script from your R environment or command line:
    ```bash
    Rscript codes/WBT_DC_pipeline.R
    ```
    This will generate results in the `results/` directory.

## License
This project is licensed under the MIT License.


