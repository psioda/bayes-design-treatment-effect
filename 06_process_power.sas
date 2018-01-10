
/*
Program Name: 06_process_power.sas
Program Creator: Matthew Psioda;

This SAS program processes the results from the second stage of the simulation 
process to identify the optimal design for the set of null and alternative 
sampling priors considered.To do this, LOESS methods are used to smooth 
the results and to interpolate the ideal values for the target number of events
and for the value of a0.

The optimal design parameters are written to the location 
&root./data/simulation_results/calculate_power (RESULTS SAS library);  The dataset created
is named results.calculated_power. The dataset contains the following variables:

[ 1] category   – description of null sampling prior used for simulation setting;
[ 2] alpha_targ – targeted Bayesian type I error rate;
[ 3] nsp        - name of the dataset that stores the null sampling prior results for the given simulation setting;
[ 4] asp        - name of the dataset that stores the alternative sampling prior results for the given simulation setting;
[ 5] shareBLH   – indicator variable for whether or not all parameters are shared between model (must be 
                  zero in current software);
[ 6] nuTarget   – desired number of events in the simulated new trials;
[ 7] a0         – value of a0 to use for power prior;
[ 8] power      - estimated bayesian power;

*/

ods html close;
proc datasets library=work noprint kill; run; quit;

/*%let root           = C:\Users\psioda\Documents\Research Papers Temp\bayesDesignTreatment\sas_programs_github;*/
/*%let resPath        = &root.\data\simulation_results\calculate_power;*/

%let root           = /nas/longleaf/home/psioda/stat-projects/bayesDesignTreatment;
%let resPath        = &root./data/simulation_results/calculate_power;

libname _all_ clear;
libname results "&resPath.";

** read in results from the first stage simulations;
data results.combined_calc_power;
 set results.calc_: indsname=ds;
  dsName = ds;
run;proc sort; by category alpha_targ nsp asp shareblh a0 nuTarget; run;

** output extra values of sample size for interpolation;
data temp;
 set results.combined_calc_power;
 by category alpha_targ nsp asp shareblh a0 nuTarget;
 output;

   nuTarget + 5;
   rejRate = .;
   a0 = .;
   output;

run;

** use PROC LOESS to smooth estimated a0 values;
proc loess data = temp;
 by category alpha_targ nsp asp shareblh ;
 model a0 = nuTarget;
 output out = smoothed_a0(drop=SmoothingParameter depvar obs) p=smooth_a0;
run;

** use PROC LOESS to smooth estimated rejection rates;
proc loess data = smoothed_a0;
 by category alpha_targ nsp asp shareblh ;
 model rejRate = nuTarget;
 output out = smoothed_rejRate(drop=SmoothingParameter depvar obs where=(rejRate=.)) p=smooth_rejRate;
run;



** calculate the optimal value of a0 for each simulation setting;
data selected;
 set smoothed_rejRate;
 by category alpha_targ nsp asp shareblh ;
 *where find(asp,'wrong','i')=0;

  power_targ = 0.80;


  diff = abs(power_targ-smooth_rejRate);
  below = (smooth_rejRate<power_targ);

  drop a0;
  rename smooth_a0 = a0;
run;

proc sort data = selected nodupkey dupout=dup;
 by category alpha_targ nsp asp shareblh diff below;
run;

proc sort data = selected out = selected_2(keep = category alpha_targ nsp asp shareblh a0 nuTarget smooth_rejRate) nodupkey;
 by category alpha_targ nsp asp shareblh ;
run;


data results.calculated_power;
 retain  category alpha_targ nsp asp shareblh ;
 set selected_2;
 rename smooth_rejRate = power;
run;

/* 

ods html;
proc contents data =  results.calculated_power varnum; run; quit; 

*/
