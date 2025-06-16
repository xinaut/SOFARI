# Yeast eQTL Data Analysis

This folder contains scripts and data files for processing the **yeast eQTL dataset** and performing statistical inference.

- **`yeast.rda`**: The raw yeast eQTL dataset.
- **`data_clean.R`**: Preprocesses `yeast.rda` and saves the cleaned data as `yeast_preprocess_data.RData`.
- **`inference_yeast.R`**: Performs the SOFARI procedure on the yeast data and produces all the **inference results** (including **Table 8**) in Appendix D of the Supplementary Material.
- **`split_yeast.R`**: Splits the cleaned data and runs both SOFARI and SOFAR methods to obtain **prediction results** in Appendix D of the Supplementary Material.
- **`yeast_inference.RData`**: Stores the inference results for all samples, allowing for direct loading without rerunning the pipeline.