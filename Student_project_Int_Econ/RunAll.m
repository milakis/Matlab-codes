% Simple file to run the programs for the paper.

clear; 
clc; 
close all;

addpath('/usr/lib/dynare/matlab') ;

%prepare Cyprus data and add it to to the original data
prepare_data_Cyprus('data_Cyprus.xlsx','Data_rebased.xlsx') ;

%solve and simulate for the structural model
Prog_Sim_Structural_NKM('pigs') ; 

%run counterfactuals
Prog_Sim_Counterfactual('pigs','cf_fiscal','fiscal counterfactual') ;
Prog_Sim_Counterfactual('pigs','cf_mp','macroprudential counterfactual') ;
Prog_Sim_Counterfactual('pigs','cf_rho','no segmentation counterfactual') ;

%produce goodness of fit table
goodness_fit_table() ;

