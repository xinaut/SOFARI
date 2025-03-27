The following scripts are used for processing the real datasets and performing inference:

- **`data_clean.R`**
   Preprocesses the raw dataset **`yeast.rda`** and outputs the cleaned data as **`yeast_preprocess_data.RData`**.
- init_code.R
- **`inference on eQTL.R`**
   Implements the **SOFARI procedure** to perform inference and generate the corresponding results.
- **`split_yeast.R`**
   Splits the yeast dataset and runs both **SOFARI** and **SOFAR** methods to obtain prediction results.
- load yeast2.RData to quickly obtain the inference results for all samples

## Yeast c Data Analysis

The following scripts are used for processing the **yeast eQTL dataset** and performing inference:

- **`data_clean.R`**
   Preprocesses the raw dataset **`yeast.rda`** and saves the cleaned data as **`yeast_preprocess_data.RData`**.
- **`init_code.R`**
   Contains initial setup code and helper functions required by the main analysis scripts.
- **`inference on eQTL.R`**
   Performs the **SOFARI procedure** to conduct inference and generate the corresponding results.
- **`split_yeast.R`**
   Splits the yeast dataset and runs both **SOFARI** and **SOFAR** methods to produce prediction results.
- **`yeast2.RData`**
   Contains the **inference results for all samples**, enabling a quick start without rerunning the full inference pipeline.