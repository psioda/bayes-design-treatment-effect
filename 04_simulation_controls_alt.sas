/*
Program Name: 04_simulation_controls_alt.sas
Program Creator: Matthew Psioda;

This SAS program generates the second stage simulation settings that were explored in the SDD example 
based on the e1684 dataset that is presented in the paper. Using the stage one results that are stored
in the folder &root./data/simulation_results/calibrate_a0 (RESULTS SAS library), this program 
first identifies the optimal value of a0 for each sample size considered in stage one subject to the 
specified Bayesian type I error control requirements. To do this, LOESS methods are used to smooth 
the results and to interpolate the ideal value of a0 for each target number of events. A dataset 
containing these results is written to the folder &root./data/simulation_results/calibrate_a0 
(RESULTS SAS library) and is named results.selected_a0. The variables included in this dataset
are named nsp, shareBLH, nuTarget, T1Erate, alpha_targ, and a0 (see below for definitions).

Using these results, the program then creates a dataset of simulation settings to be explored in 
the second stage of the SSD procedure (named controls.simulation_controls_alt). The SAS dataset of 
simulation settings is written to the location: &root./data/simulation_controls (CONTROLS SAS 
library). The dataset contains the following variables:
[ 1] category       – description of null sampling prior used for simulation setting;
[ 2] nsp            - name of the dataset that stores the null sampling prior results for the given simulation setting;
[ 3] asp            - name of the dataset that stores the alternative sampling prior results for the given simulation setting;
[ 4] shareBLH       – indicator variable for whether or not all parameters are shared between model (must be 
                      zero in current software);
[ 5] nuTarget       – desired number of events in the simulated new trials;
[ 6] T1Erate        – estimated Bayesian type I error rate associated with number of events and null sampling prior
[ 7] alpha_targg    – targeted Bayesian type I error rate
[ 8] a0             – value of a0 to use for power prior
[ 9] node_idx       – Integer value that identifies the compute node to be used for a given set of simulations;
[10] nPerLoop       – number of datasets to simulate in a given loop cycle – this number cannot be too large due to 
                      system resources. A value of 500 is used in the current setup. The more datasets that are simulated at 
                      once, the more memory or disk space needs to be available to SAS to process them;
[11] inner_idx      – observation counter for simulation settings to be performed by a given compute node;
[12] seed           – random number seed for new trial data generation;
[13] parmSampleSeed – random number seed for drawing samples from discrete approximation to the 
                      sampling priors;
*/


proc datasets library=work noprint kill; run; quit;
ods html close;
ods listing close;

%let root        = /nas/longleaf/home/psioda/stat-projects/bayesDesignTreatment;
%let contPath    = &root./data/simulation_controls;
%let resPath     = &root./data/simulation_results/calibrate_a0;

/*%let root        = C:\Users\psioda\Documents\Research Papers Temp\bayesDesignTreatment\sas_programs_github;*/
/*%let contPath    = &root.\data\simulation_controls;*/
/*%let resPath     = &root.\data\simulation_results\calibrate_a0;*/


libname _all_ clear;
libname controls "&contPath.";
libname results "&resPath.";


%let number_of_nodes = 600;
%let nPerLoop        = 500;
%let nLoops          = 200;

** read in results from the first stage simulations;
data results.combined_calibrate_a0;
 set results.calibrate:;
  by sp shareblh nuTarget a0;
  category = catx("-",strip(nuTarget),strip(shareblh));
run;proc sort; by sp nuTarget shareblh; run;

** output extra values of a0 for interpolation;
data temp;
 set results.combined_calibrate_a0;
 by sp nuTarget shareblh;

  output;

  if last.shareblh then do a0 = 0.0 to 1.0 by 0.005;
   rejRate = .;
   output;
  end;
run;


** use PROC LOESS to smooth estimated rejection rates;
proc loess data = temp;
 by sp nuTarget shareblh;
 model rejRate = a0;
 output out = smoothed_a0(where=(rejRate = .)) p=smooth_rejRate;
run;


** calculate the optimal value of a0 for each simulation setting;
data selected;
 set smoothed_a0;
 by sp nuTarget shareblh;

  alpha_targ = 0.025;

  diff  = abs(alpha_targ-smooth_rejRate);
  above = (smooth_rejRate>alpha_targ);
  if sp in ("sPrior.default_null" "sPrior.general_null" "sPrior.freq_null") then output;
run;

proc sort data = selected nodupkey dupout=dup;
 by sp shareblh nuTarget alpha_targ diff above;
run;

proc sort data = selected out = selected_2(keep = sp nuTarget shareblh alpha_targ a0 smooth_rejRate) nodupkey;
 by sp shareblh nuTarget alpha_targ;
run;

proc sort data = selected_2;
 by sp shareblh alpha_targ nuTarget;
run;

** use PROC LOESS to smooth estimated a0 values;
proc loess data = selected_2;
 by sp shareblh alpha_targ;
 model a0 = nuTarget;
 output out = selected_2b(drop=a0 rename=(smooth_a0=a0)) p=smooth_a0;
run;

data results.selected_a0;
 set selected_2b;
 drop SmoothingParameter DepVar Obs;
 rename sp = nsp smooth_rejRate = T1Erate;
run;

data selected_3;
 set selected_2b(rename=(sp=nsp));
 length asp $50.;
 
 asp = 'sPrior.default_alt'; output;
 asp = 'sPrior.general_alt'; output;
 asp = 'sPrior.point_mass_alt'; output;
 asp = 'sPrior.freq_null'; output;

run; proc sort; by asp nsp shareblh alpha_targ nuTarget; run;
 
data selected_4;
 retain category nsp asp;
 set selected_3;

  length category $50.;
       if nsp = 'sPrior.default_null' then Category = 'Default Null Sampling Prior';
  else if nsp = 'sPrior.general_null' then Category = 'General Null Sampling Prior';
  else if nsp = 'sPrior.freq_null'    then Category = 'Frequentist-Like Null Sampling Prior';

  drop SmoothingParameter DepVar Obs;
run;


data controls.simulation_controls_alt;
 set selected_4;

  nPerLoop = &nPerLoop.;


  retain node_idx 0;
  node_idx + 1;
  if node_idx > &number_of_nodes. then node_idx = 1;

 do inner_idx = 1 to &nLoops.;
   seed           = round(1 + rand('uniform')*2**30);
   parmSampleSeed = round(1 + rand('uniform')*2**30);
    output;
 end;

 rename smooth_rejRate = T1Erate;
run;


proc sort data = controls.simulation_controls_alt; 
 by category nsp asp shareblh nuTarget alpha_targ a0 ; 
run;

/*
ods html;
proc contents data = controls.simulation_controls_alt varnum; run;
*/
