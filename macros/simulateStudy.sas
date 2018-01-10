/*
 This SAS macro is used simulate data for the new trial. The macro simulates data from a proportional hazards model with piecewise constant
 baseline hazard.
 The macro requires the following parameters:
  
  MACRO VARIABLES REQUIRED:
   [1] trueHazComp     = list of the number of baseline hazard components in each stratum - each stratum separated by a pipe;
   [2] dsSamples       = dataset that contains the samples for the chosen sampling prior;
   [3] parmSampleSeed  = seed used by PROC SURVEYSELECT when sampling parameter values from sampling prior dataset using 
                         unrestricted random sampling;
   [4] numSimulations  = number of datasets to be simulated;
   [5] seed            = seed value used in random number generator when simulating new trial data;
   [6] dsBreaks        = dataset that stores the baseline hazard change points - required as input;
   [7] regParmNumber   = number of hazard ratio regression parameters;
   [8] stratLevels     = number of strata;
   [9] regParmDist     = distribution for each covariate being simulated for the hazard ratio regression model - distributions 
                         should be valid calls to the RAND function and should be separated by a pipe;
   [10] stratProbs     = probabilities of subject allocation to strata - list all but last probability as it will be 
                         computed as 1 - sum(p_s) - probabilities should be separated by commas;
   [11] enrollDist     = distribution for enrollment - should be valid call to RAND function or a constant;
   [12] censorDist     = distribution for censorship - should be valid call to RAND function or a constant;
   [13] nuTarget       = targeted number of events for new trial;
   [14] dsOut          = dataset to store simulated new trial data;
*/




%macro simulateStudy;
%let maxComp = 0;
%do s = 1 %to &stratLevels.;
  %if %eval(&maxComp < %scan(&trueHazComp.,&s,|)) %then %let maxcomp = %scan(&trueHazComp.,&s,|);
%end;


proc surveyselect data = &dsSamples. out = __post_samples_selected__(drop=NumberHits) seed=&parmSampleSeed.
          outhits method=urs n=&numSimulations. noprint; run; quit;

data __temp__;
 call streaminit(&seed.);
 set __post_samples_selected__;
  if _n_ = 1 then set &dsBreaks.;

 simStudy = _n_;

 array x(&regParmNumber.);
 array b(&regParmNumber.) %do p = 1 %to &regParmNumber.; b&p. %end; ;


 array h(&stratLevels.,&maxComp.) %do s = 1 %to &stratLevels.; %do c = 1 %to &maxComp.; h&s.&c. %end; %end; ;
 array t(&stratLevels.,&maxComp.) %do s = 1 %to &stratLevels.; %do c = 1 %to &maxComp.; t&s.&c. %end; %end;;


 do subID = 1 to &numSubjects.;
  ** calculate subject hazard ratio regression function;
  logHR = 0;
  %do j = 1 %to &regParmNumber.;
   x(&j.) = %scan(&regParmDist.,&j,|);;
   logHR = logHR + x(&j.)*b(&j);
  %end;
  HR = exp(logHR);

  *** simulate stratum for subject;
  riskStratum = rand('table',&stratProbs.);

  ** simulate enrollment time;
  enrollTime = round(&enrollDist.,0.0001);

  ** simulate censorship time;
  censTime = round(&censorDist.,0.0001);

  ** Simulate data from Piecewise-exponential model;
  k       = 1;
  obsTime = 0;
  stop    = 0;
  do while(stop=0);
   haz = h(riskStratum,k)*HR;
   obsTime   = obsTime + rand('exponential')/haz;

   if obsTime > t(riskStratum,k) and t(riskStratum,k) > . then do; obsTime = t(riskStratum,k); k+1; end;
   else stop = 1;
  end;
  obsTime = round(obsTime,0.0001);
  if obsTime > censTime then do; obsTime = censTime; censor = 1; end;
  else censor = 0;

  elapsedTime = enrollTime + obsTime;

  output;
 end;
 
 keep obsTime enrollTime elapsedTime x: censor riskStratum subID simstudy;
run;

proc sort data = __temp__ out = __events__;
 by simstudy elapsedTime;
 where censor = 0;
run;

data __events1__;
 set __events__ ;
 by simstudy;
 
  if first.simstudy then eventCounter = 0; eventCounter +1;
  
  if eventCounter = &nuTarget. or (last.simstudy and eventCounter < &nuTarget.) then output;
run;

data __events2__;
 set __events1__;
  targetTime = elapsedTime + 0.00005;
  keep simstudy targetTime;
run;

data &dsOut.;
 merge __temp__ __events2__;
 by simstudy;
 
 if enrollTime >= targetTime then delete;
 else if elapsedTime > targetTime then do;
   obsTime = targetTime - enrollTime;
   censor = 1;
   elapsedTime = enrollTime + obsTime;
 end;

run;


 proc datasets library=work noprint;
  delete __temp__ __events__ __events1__ __events2__;
 run;
 quit; 
%mend simulateStudy;
