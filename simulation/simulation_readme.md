## Reproducing the Simulation Studies

To reproduce the simulation studies presented in the paper, please refer to the following files:

- **`result_sum.R`**  
  This is the **main simulation script**, which performs **parallel computing** and **outputs the simulation results**.

- **`main of secure and sofar.R`**  
  This file contains the main functions for the **three simulation examples** in the paper:
  - `pare_secure`: **Simulation Example 1**
  - `pare_sofar`: **Simulation Example 2**
  - `pare_sofar_nonsparse`: **Simulation Example 3**

- **`functions.R`**  
  Contains a collection of **basic utility functions** used across simulation scripts.

- **`func_nearly.R`**  
  Provides functions to **compute the distribution of strongly orthogonal factors**.

- **`func_weakly.R`**  
  Provides functions to **compute the distribution of weakly orthogonal factors**.

- **`sim.R`**  
  Estimates the **covariance matrix of random errors**  using **adaptive thresholding**.

- **`STRS.R`**  
  Estimates the **rank of multivariate response regression models**.

- **`helpers.nodewise.R`**  
  Estimates the **precision matrix** (i.e., the inverse of the covariance matrix).

