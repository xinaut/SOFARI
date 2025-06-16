# Economic Data Analysis

This folder contains scripts and data files for analyzing the **Federal Reserve Economic Data** used in the real data section of the paper.

- **`2015-04.csv`**: The raw economic dataset.
- **`economic_data_process.R`**: Preprocesses `2015-04.csv` and outputs the cleaned covariate and response matrices `final_XX.csv` and `final_YY.csv`.
- **`inference_economic.R`**: Performs the inference procedure using the preprocessed data. It produces the result files `final_Uk.csv`, `final_theta.csv`, `final_Sigmae.csv`, `dU.csv`, and `V.csv`. It also generates the results for **Table 2** and **Table 3** in the paper.
- **`Table 7.R`**: Generates the results for **Table 7** in the paper.
- **`plot_economic.R`**: Visualizes the features selected by the inference procedure and produces **Figure 2**.
- **`realdata1.RData`**: Stores the preprocessed data and full inference results, allowing for a quick start without rerunning the full pipeline.
- **`realdata2.RData`**: Contains rearranged data and results used specifically for generating Figure 2.