## Economic Data Analysis

The following scripts and data files are used for processing and analyzing the **Federal Reserve Economic Data**:

- **`real_data_process.R`**
   Preprocesses the raw data in **`2015-04.csv`** and outputs the datasets **`final_XX.csv`** and **`final_YY.csv`**.
- **`infer of real data.R`**
   Performs the **inference procedure** and generates the corresponding inference results.
- **`plot_realdata.R`**
   Visualizes the **features selected by the inference procedure**.
- **`realdata1.RData`**
   Contains the **preprocessed data** and the **inference results**, allowing for a quick start without rerunning the full pipeline.
- **`realdata2.RData`**
   Contains the **inference results** and the **rearranged data** used specifically for plotting.