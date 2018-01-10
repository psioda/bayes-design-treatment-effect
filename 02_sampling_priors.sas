/*
Program Name: 02_sampling_priors.sas
Program Creator: Matthew Psioda;

This SAS program can be run to construct discrete approximations to the sampling priors used in the paper. 

The program first uses the %COLLAPSE macro to collapse the subject-level e1684 dataset to a dataset of 
sufficient statistics. The model fitting SAS macros (%FIT_MCMC and %FIT_APPROX) expect this format. This 
data reduction allows the model to be fit using PROC GENMOD using the connection between PH regression 
and Poisson regression.

Using the reduced historical dataset, the %FIT_MCMC macro is repeatedly called and the results are 
post-processed. This sas program includes the macro definition for the %EST_SAMP_PRIOR macro. The 
%EST_SAMP_PRIOR is nothing more than a wrapper to the %FIT_MCMC macro that does some post-processing 
of the results.

This program generates the discrete approximations for the default, truncated, and point mass 
sampling priors that were used in the paper. The macro calls that generate the discrete approximations 
begin just after the macro definition for the %EST_SAMP_PRIOR.

This program reads in the raw data from the folder &root./data/raw_data (RAW SAS library) 
and writes out the sampling prior datasets to the folder &root./data/sampling_priors (RES SAS library);

The program reads in the %COLLAPSE and %FIT_MCMC macros from the folder &root./macros;

*/


proc datasets library=work noprint kill; run; quit;


%let root        = /nas/longleaf/home/psioda/stat-projects/bayesDesignTreatment;
%let rawDataPath = &root./data/raw_data;
%let resDataPath = &root./data/sampling_priors;
%let macroPath   = &root./macros;

%include "&macroPath./collapse.sas";
%include "&macroPath./fit_mcmc.sas";

libname _all_ clear;
libname raw "&rawDataPath.";
libname res "&resDatapath.";

data stage4;
 set raw.e1684_stage4;
run;


%let dsIn         = stage4;                   ** name of dataset to collapse to sufficient statistics;
%let dsOut        = raw.stage4_condensed;     ** name of dataset to write out;

%let writeBreaks  = 1;                        ** indicator of whether or not dataset should be created to store hazard change points;
%let dsBreaks     = raw.hazard_breaks;        ** name of dataset that will contain hazard break points;


%let stratLevels  = 2;                         ** number of strata;
%let fitHazComp   = 4 | 3;                     ** number of hazard components in each stratum;
%let fitHazBreaks = 0.15342 0.46301 2.29589 |  
                    0.59726 2.87397;           ** hazard change points;

%collapse; 





/* empty historical dataste for use with the %FIT_MCMC macro; */
data hist; set &dsOut.(obs=0); run;





/* This macro post-processes results for the %FIT_MCMC macro
   The %FIT_MCMC macro is called within this post-processing macro
*/
%macro est_samp_prior;
%fit_mcmc;

proc sort data = &dsOutParmEst._summary sortseq=linguistic(numeric_collation=on) out = res.&dsOutParmEst._summary;
 by parameter;
run;

data &dsOutParmEst.;
 set &dsOutParmEst.(drop=HR postProb log: iteration simStudy iteration);
  sort = rand('uniform');
run; proc sort data=&dsOutParmEst. out=&dsOutParmEst.(drop=sort); by sort; run;



** count number of strata;
%let numStrat = %eval(%sysfunc(count(&fithazComp,|))+1);

** count maximum number of compoents;
%let maxComp = 0;
%do s = 1 %to &numStrat.;
  %if %eval(&maxComp < %scan(&fitHazComp.,&s,|)) %then %let maxcomp = %scan(&fitHazComp.,&s,|);
%end;
%put maxcomp = &maxcomp.;

data temp;
  %let var = 1;
 %do s = 1 %to &numStrat.;
 %do i = 1 %to &maxComp.; 
  h&s.&i.  = .;
  %let var = %eval(&var + 1);
 %end;
 %end;

data &dsOutParmEst.;
 set &dsOutParmEst.;

 rename gamma = b1

        %let var = 1;
        %do s = 1 %to &numStrat.;
        %do i = 1 %to  %scan(&fitHazComp.,&s,|); 
          lambda&var. = h&s.&i. 
          %let var = %eval(&var + 1);
        %end; 
        %end;;
 run;

 data res.&dsOutParmEst.;
  set temp(in=a) &dsOutParmEst.;
   if not a;
 run;

%mend est_samp_prior;


** initial priors;
%let gammaPrior   = normal(0,var=1e5);
%let lambdaPrior  = gamma(1e-5,iscale=1e-5);

** new and historical sets to analyze;
%let dsSim        = raw.stage4_condensed;
%let dsHist       = hist;

%let shareBLH     = 1;      ** indicator of whether nuisance parameters are shared;
%let a0           = 0.00;   ** a_0 value for partial-borrowing power prior;
%let dsOutDIC     = DIC;    ** output dataset to store model DIC estimate;


%let seedValue    = 325211;           ** seed for PROC MCMC;
%let dsOutParmEst = e1684_posterior;  ** name of dataset of posterior samples;
%let mcmcOptions  = nmc=500000 nbi=500 nthin=5 seed=&seedValue. postout=&dsOutParmEst.; ** PROC MCMC options;
%est_samp_prior;

proc means data = res.E1684_posterior noprint nway;
 var h: b:;
 output out = res.point_mass_alt(drop=_:) mean= ;
run;

%let seedValue    = 624562;
%let dsOutParmEst = default_null;
%let gammaPrior   = normal(0,var=1e5,lower=0);
%let mcmcOptions  = nmc=500000 nbi=500 nthin=5 seed=&seedValue. postout=&dsOutParmEst.; 
%est_samp_prior;

%let seedValue    = 684686;
%let dsOutParmEst = general_null;
%let gammaPrior   = normal(0,var=1e5,lower=0,upper=0.08215);
%let mcmcOptions  = nmc=500000 nbi=500 nthin=5 seed=&seedValue. postout=&dsOutParmEst.; 
%est_samp_prior;

%let seedValue    = 8658654;
%let dsOutParmEst = default_alt;
%let gammaPrior   = normal(0,var=1e5,upper=0);
%let mcmcOptions  = nmc=500000 nbi=500 nthin=5 seed=&seedValue. postout=&dsOutParmEst.; 
%est_samp_prior;

%let seedValue    = 7646442;
%let dsOutParmEst = general_alt;
%let gammaPrior   = normal(0,var=1e5,lower=-0.48997,upper=-0.04107);
%let mcmcOptions  = nmc=500000 nbi=500 nthin=5 seed=&seedValue. postout=&dsOutParmEst.; 
%est_samp_prior;

%let seedValue    = 452426;
%let dsOutParmEst = freq_null;
%let gammaPrior   = normal(0,var=1e-22);
%let mcmcOptions  = nmc=500000 nbi=500 nthin=5 seed=&seedValue. postout=&dsOutParmEst.; 
%est_samp_prior;

