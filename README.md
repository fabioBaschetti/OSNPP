# Replication Package

**An Optimal Energy Production Problem with Energy Source Switching and Load Following Nuclear Power Plants**

## Software Requirements

- **Software:** MATLAB
- **Version:** R2023b
- **Toolboxes:**
  - Parallel Computing Toolbox
  - Statistics and Machine Learning Toolbox
  - Econometrics Toolbox
- **Operating System:** Windows*

## Folder Structure

```text
.
├── data_ITA/                         store data for ...
│   ├── renewables/                   ... renewable production 
│   └── total_load/                   ... total demand
│                                     process that data
│                                     and performs ML estimation of an AR (/ OU) process for the residual demand
│
├── OSNPP_clsdmkt/                    experiments in a closed market
│
├── OSNPP_openmkt_pricemaker/         experiments an an open  market where the local producer is a price maker
│   └── policy_table/                 store policy table for the current parameters (no need to run the solver^)
│
├── OSNPP_openmkt_pricetaker/         experiments an an open  market where the local producer is a price taker
│   └── policy_table/                 store policy table for the current parameters (no need to run the solver^)      
│
├── paths/                            store some simulated paths for the residual demand (both local and global)
                                      as needed for producing the Figures in the paper
│
├── summary_table/                    store all scripts genrating Table 1 in the paper
│
└── utils/                            auxiliary functions for solving the optimal switching problem, simulate a
                                      controlled diffusion and plot the results, both ...
    ├── 2D/                           ... in the closed-economy case
    └── 3D/                           ... in the open  -economy case
```

## Running the Code

To generate **Table 1**:

```matlab
results = run_summary_table;
```

To run individual experiments (with figures):

```matlab
main.m
```

## Data

Data in the `data_ITA` folder are publicly available at **Terna Download Center**:

https://dati.terna.it/download-center

## Folder-Specific Instructions

### `data_ITA`

- **renewables.m**
  - Processes renewables data in `renewables/`
  - Produces `renewables.mat`

- **total_load.m**
  - Processes total load data in `total_load/`
  - Produces `total_load.mat`

- **mle_Y.m**
  - Uses `total_load.mat` and `renewables.mat`
  - Builds a residual demand time series
  - Fits an Ornstein–Uhlenbeck process
  - The estimated parameters are used by the main scripts

- **plot_IT_vs_EU.m**
  - Produces **Figure 10** in the paper
  - Here **EU** is used as short-hand notation for a fictitious system that includes Italy.

### `OSNPP_clsdmkt`

The user can manually modify the settings to reproduce the sensitivity analyses presented in the paper and change the time horizon.

## Notes

- * Some simulations were run on macOS. The corresponding outputs are stored in the `paths` folder and accessed by the `main.m` scripts.*
- ^ Running the **3D solver** over a **7-day horizon** takes approximately **35 minutes** on a standard laptop with **8 CPU cores**.
- ^ Running the **2D solver** is substantially faster and can be done online (this is why we do not store the optimal policy for that).
