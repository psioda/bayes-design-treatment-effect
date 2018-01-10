/*
Program Name: 03_calibrate_a0.sas
Program Creator: Matthew Psioda;

This SAS program is designed to complete stage one of the sample size determination simulations. That is, 
this program is designed to compute the Bayesian type I error rate for given combinations of the number 
of events targeted in the new trial and values of a0. The program simulates new trial datasets based on 
input design settings, analyzes the data using the partial-borrowing power prior, and then computes 
estimates of the Bayesian type I error rate using the results from the simulated trials. 

The program first uses the following macros:
[1] %SIMULATE_STUDY --> Simulates data for the new trial;
[2] %COLLAPSE       --> Reduces the new trial data to the sufficient statistics;
[3] %FIT_APPROX     --> Analyzes the new trial data using the power prior;

The macros are expected to reside in the location &root./macros;

This program reads are writes data to the following locations:
[1] &root./data/raw_data --> location from which raw data is read (RAW SAS library);
[2] &root./data/sampling_priors --> location for sampling prior datasets (SPRIOR SAS library);
[3] &root./data/simulation_controls --> location from which simulation settings are read (CONTROLS SAS library);
[4] &root./data/simulation_results/calibrate_a0 --> location to which simulation analysis results are written (RESULTS SAS library);

This program is designed to be run multiple times in parallel (likely using a command line interface) 
where the users passes the SAS program a macro variable value through the automatic SYSPARM macro variable. 
The value of the SYSPARM macro variable identifies the compute node for a computing cluster environment. 
For each value of the SYSPARM macro variable, a different set of simulation settings are used based on 
the controls.simulation_controls_null dataset created by a previously run program.

This program includes the macro definition for the %LOOP macro that reads in the simulation setting 
observations from the controls.simulation_controls_null dataset that are assigned to the compute node 
in question and then loops over the simulation settings to perform the specified number of simulations 
for each simulation setting (i.e., observation in the controls.simulation_controls_null dataset).

*/ 


%let root           = /nas/longleaf/home/psioda/stat-projects/bayesDesignTreatment;
%let rawDataPath    = &root./data/raw_data;
%let macroPath      = &root./macros;
%let sPriorDataPath = &root./data/sampling_priors;
%let contPath       = &root./data/simulation_controls;
%let resPath        = &root./data/simulation_results/calibrate_a0;

%include "&macroPath./simulateStudy.sas";
%include "&macroPath./collapse.sas";
%include "&macroPath./fit_approx.sas";

libname _all_ clear;
libname raw "&rawDataPath."       access=read;
libname sPrior "&sPriorDataPath." access=read;
libname controls "&contPath."     access=read;
libname results "&resPath.";


proc datasets library=work kill noprint;quit;

%macro loop;
%if &debug=1 %then %let node_id=1;
data time_track;
 start_time = time();
run;

data _null_;
 call symput('dsLabel',put(&node_id.,z5.));
run;

%** read in design parameters for simulation;
data work.__controls__;
 set controls.simulation_controls_null;
 where node_idx = &node_id.;

 if &debug = 1 then do;
  if _n_> 5 then delete;
  nPerLoop = 100;
 end;

 row_idx = _n_;
 call symput('nSimSettings',strip(put(_n_,best.)));
run;

%** read in baseline hazard change points;
data _null_;
 set &dsBreaks.;
 array t(*) t:;

 length fitHazBreaks $5000.;
 fitHazBreaks = "";
  do j = 1 to dim(t);
    if t(j) >. then do; fitHazBreaks = catx(" ",fitHazBreaks,strip(put(t(j),best.))); end;
    else if t(j) = . and substr(strip(fitHazBreaks),length(strip(fitHazBreaks)),1) ^= '|' then do;
     fitHazBreaks = catx(" ",fitHazBreaks,"|");
    end;
  end;

   if substr(fitHazBreaks,length(fitHazBreaks),1) = '|' then substr(fitHazBreaks,length(fitHazBreaks),1) = ' ';
   call symput('fitHazBreaks',strip(fitHazBreaks));
 run;

%** loop over different design parameters settings for simulations;
%do loop     = 1 %to &nSimSettings.;
%if %sysfunc(mod(&loop,10))=1 %then %put loop &loop of &nSimSettings.;

  %** create macro variables needed for %SIMULATE_STUDY macro and %FIT_APPROX macro;
  data _null_;
   set work.__controls__(drop=inner_idx);
   where row_idx = &loop.;
   call symput('a0',strip(put(a0,best.)));

   call symput('nuTarget',strip(put(nuTarget,best.)));
   call symput('numSubjects',strip(put(nuTarget*3,best.)));

   call symput('shareBLH',strip(put(shareBLH,best.)));

   call symput('dsSamples',strip(sp));

   call symput('seed',strip(put(seed,30.)));
   call symput('parmSampleSeed',strip(put(parmSampleSeed,30.)));

   call symput('numSimulations',strip(put(nPerLoop,best.)));

  run;

  ** data set name for simulated data;
  %let dsOut     = work.simData;

  ** simulate data;
  %let trueHazComp = &fitHazComp_sim.;
  %simulateStudy;

  ** collapse data;
  %let dsIn         = work.simData;
  %let dsOut        = work.simDataCondensed;

  %let writeBreaks  = 0;
  %if &shareBLH ^= 0 %then %abort; %* functionality for SHAREBLH=1 removed;
  %let fitHazComp = &fitHazComp_sim.;
  %collapse;


  %* fit the model using weighted maximum likelihood approximation with %FIT_APPROX macro;
  %let dsSim         = work.simDataCondensed;
  %let dsHist        = raw.stage4_condensed;
  %let dsOutParmEst  = simParmEst_Approx;
  %fit_Approx;

  %* extract needed results from analysis output;
  data work.__Results;
   set &dsOutParmEst. end=last;
   where Parameter = 'POSTPROB';

   row_idx = &loop.;
  run;

  %* append results in ongoing loop;
  proc append data = work.__Results base=work.Results force; run; quit;

  %* scour work library of results from loop iteration;
  proc datasets library=work noprint;
   delete __Results;
  quit;

%end;

data work.Results2;
 merge work.__controls__(in=a) work.Results(in=b) ;
 by row_idx;
  if a and b;
run;

%* write out results to permanent location;
data results.calibrate_a0_&dsLabel.;
 set work.Results2;
 by sp shareBLH nuTarget a0;

  if first.a0 then do; rejRate=0; N=0; end;
  
  retain rejRate Den 0;
  if estimate > 0.975 then rejRate + 1;
                           N + 1;

  if last.a0 then do;
    rejRate = rejRate / N;
    output;
  end;

  keep sp shareBLH nuTarget a0 N rejRate;
run;


data time_track;
 set time_track;
 stop_time = time();

 elapsed_time = (stop_time-start_time)/60;

 put start_time= stop_time= elapsed_time=;
run;



%mend;


** dataset names baseline hazard breaks;
%let dsBreaks  = raw.hazard_breaks;

** stratification variables;
%let stratLevels   = 2; 
%let stratProbs    = 0.50;

** number of baseline hazard components in fitted model / simulated data;
%let fitHazComp_sim = 4 | 3;

** number of regression parameters and covariate distribution;
%let regParmNumber = 1;
%let regParmDist   = rand('Bernoulli',0.5);

** enrollment and censorship Distributions;
%let enrollDist    = rand('uniform')*4;
%let censorDist    = 100;

%*in general sysparm is designed to be set from the command line so that this program
  can be used for an array computing job on a high performance computing cluster. Users
  can always set sysparm themselves or directly set node_id;

%let node_id = %scan("&sysparm", 1, " ");
%put &=node_id;

** indicator for whether to use debug mode. If debug mode is on then
   at most five design parameter settigs will be evaluated and 
   nPerLoop will be set to 100;
%let debug=0;
option nonotes;
%loop;
option notes;

