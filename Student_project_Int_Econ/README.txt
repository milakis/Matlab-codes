Instructions to run the code for the term paper "Martin and Philippon (2017) model: framework to analyze the eurozone crisis in Cyprus?" by Liudmila Kiseleva

- RunAll.m: runs all the programs for the paper;

- prepare_data_Cyprus.m:
  - uploads the row data for Cyprus from the file data_Cyprus.xlsx;
  - rescales it as described in the calibration part;
  - adds new series to the original file Data_rebased.xlsx to obtain the final data data.xlsx;

- Prog_Sim_Structural_NKM.m: solves and produces simulations for the structural model;

- polf2ormap.m: maps time t-1 determined state variables into time t variables, given second order policy functions from dynare;

- sim_bg.m: simulates a series for government debt;

- figures.m: produces the figures;

- SOEfixgovSB_str.mod: dynare file called to solve the model;

- Prog_Sim_Counterfactual.m: produces counterfactual simulations;

- figures_cf.m: produces figures for counter exercises;

- goodness_fit_table.m: produces goodness of fit measures of benchmark structural model and saves them to the file fit_table_struct_demean.csv

* I use Matlab version R2019b and Dynare version 4.2.4



