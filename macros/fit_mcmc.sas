/*
 This SAS macro is used fit the PH model using the power prior via MCMC. 
 The macro requires the following parameters:
  
  MACRO VARIABLES REQUIRED:
   [1] dsSim           = dataset of simulated new trials after reduction to sufficient statistics; This dataset needs to have been 
                         created using the %COLLAPSE macro or needs to be in the same structure as the dataset produced by that macro.
   [2] dsHist          = historical dataset after reduction to sufficient statistics; This dataset needs to have been 
                         created using the %COLLAPSE macro or needs to be in the same structure as the dataset produced by that macro.
                         If there is no historical dataset, this dataset should simply have zero observations.
   [3] shareBLH        = indicator variable for whether or not baseline hazard is shared between studies; The value should be set to
                         a value of 0 for the partial-borrowing power prior.
   [4] a0              = value of a0 for the power prior;
   [5] dsOutParmEst    = prefix to use for datasets storing results from PROC MCMC;
   [6] mcmcOptions     = PROC statement options for PROC MCMC; This option can be used to set the MCMC seed;
   [7] lambdaPrior     = prior for baseline hazard parameters using PROC MCMC syntax;
   [8] gammaPrior      = prior for the treatment effect using PROC MCMC syntax;
*/

%macro fit_mcmc;

 proc sql noprint;
  select max(simStudy) into :numSimStudies
  from &dsSim.;
 quit;

proc sql noprint;
  select max(riskStratum) into :numRiskStratum
  from &dsSim.;
 quit;
 

data __tempHist__;
 set &dsHist.(in=a);

 if &shareBLH = 0 then riskStratum = riskStratum + &numRiskStratum.;

  do simStudy = 1 to &numSimStudies.;
   output;
  end;
run; 

data __tempData__;
 set &dsSim. __tempHist__(in=a);

 if a then study = 0; else study  = 1;
 if a then wgt = &a0.; else wgt = 1.0;
run;

proc sort data = __tempData__ out = __hazCounter__(keep = riskStratum interval) nodupkey;
 by riskStratum interval;
run;

data __hazCounter__;
 set __hazCounter__;
  
   
  retain haz_id 0;
  haz_id+1;
  
  length Parameter $25.;
  PARAMETER = "LAMBDA"||strip(put(haz_id,best.));
run;

proc SQL noprint undo_policy=none;
 create table __tempData__ as 
 select a.*,b.haz_id
 from __tempData__ as a left join __hazCounter__ as b
  on a.riskStratum=b.riskStratum and a.interval=b.interval
  order by simStudy, riskStratum, interval;
quit;

proc sql noprint;
 select max(haz_id) into :numHazComp from __tempData__;
quit;
%let numHazComp = %sysfunc(compress(&numHazComp.));

sasfile __tempData__ load;
ods _all_ close;
ods output PostSummaries = &dsOutParmEst._summary DIC=DIC PostIntervals  = &dsOutParmEst._intervals;
proc mcmc data = __tempData__ plots=(none) monitor=(lambda1-lambda&numHazComp. gamma HR postProb ) statistics=ALL &mcmcOptions. dic;
 by simStudy;

 array lambda[&numHazComp.];
 array logLambda[&numHazComp.];

 beginnodata;
  do k = 1 to dim(lambda);
   logLambda[k] = log(lambda[k]);
  end;
  
   HR = exp(gamma);
   postProb = (HR<1.0);
 endnodata;

 parms lambda: / slice;
 parms gamma / slice;

 prior lambda: ~ &lambdaPrior.;
 prior gamma ~ &gammaPrior.;

 likelihood = wgt*(
                    numEvents*gamma*x1 + numEvents*logLambda[haz_id] - riskTime*lambda[haz_id]*exp(gamma*x1)
                  );
 model general(likelihood);
run;
sasfile __tempData__ close;

data &dsOutParmEst._summary;
  merge &dsOutParmEst._summary &dsOutParmEst._intervals;
  parameter = upcase(Parameter);
 run;


data &dsOutParmEst._summary;
 length Parameter $25.;
  set &dsOutParmEst._summary;
  parameter = upcase(Parameter);
 run;
 
 proc sql noprint undo_policy=none;
  create table &dsOutParmEst._summary as
  select a.*,b.riskStratum,b.interval
  from &dsOutParmEst._summary as a left join __hazCounter__ as b
   on a.parameter = b.parameter;
 quit;
 
 
 data &dsOutParmEst._summary;
  set &dsOutParmEst._summary;
    if riskStratum >. then parameter = 'LAMBDA_'||strip(put(riskStratum,best.))||"_"||strip(put(interval,best.));
    keep simStudy N mean stdDev p50 parameter HPD:;
 run; proc sort; by simStudy parameter; run;

/*
  proc transpose data = &dsOutParmEst. out = PostSummariesAll(drop=_:);
   by simSTudy;
   id Parameter;
   var mean;
  run;
*/
 proc datasets library=work noprint;
  delete __tempData__ __tempHist__ __hazCounter__ ;
 run;
 quit;  

%mend fit_mcmc;