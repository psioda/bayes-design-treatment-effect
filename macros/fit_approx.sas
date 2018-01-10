/*
 This SAS macro is used fit the PH model using the power prior via weighted maximum likelihood approximation. 
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
   [5] dsOutParmEst    = prefix to use for datasets storing results from PROC GENMOD;
*/


%macro fit_Approx;


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
run; proc sort; by simStudy; run;

data __tempData__;
 set &dsSim. __tempHist__(in=a);
 by simStudy;
 if a then wgt = &a0.; else wgt = 1.0;
 
 
 ** Fudge Factors;
  if riskTime = 0 then riskTime = 1e-7;
  if numEvents = 0 then numEvents = 1e-7;
  
  logRiskTime = log(riskTime);

 riskInt = strip(put(riskStratum,z3.))||'-'||strip(put(interval,z3.));
run;

 ods select none;
 ods output ParameterEstimates = &dsOutParmEst. ;
 proc genmod data = __tempData__;
  by simStudy;
  weight wgt;
  class riskInt(ref='001-001');
  model numEvents = riskInt x: / dist=poisson link=log offset=logRiskTime;
 run;
 quit;

 data &dsOutParmEst.;
  length parameter $25.;
  set &dsOutParmEst.;

   parameter = upcase(parameter);
  
  retain int;
  if parameter = 'INTERCEPT' then do; int = estimate; end;
  if find(parameter,'RISKINT','i') then do;
      estimate = exp(estimate + int);
	  hazCounter+1;
	  parameter = 'LAMBDA_'||strip(put(input(scan(Level1,1,'-'),best.),best.))||"_"||strip(put(input(scan(Level1,2,'-'),best.),best.));
  end;
  

  
  if parameter = 'X1' then PARAMETER = 'GAMMA';
  if parameter not in ("SCALE" "INTERCEPT") then output;
   
  if PARAMETER = 'GAMMA';
  
  est = estimate;
  std = stderr;
  
  parameter= 'HR';
  estimate = exp(est);
  stderr   = sqrt(std**2*est**2);
  output;
  
  parameter= 'POSTPROB';
  estimate  = 1-CDF('normal',est/std);
  stderr = .;
   output;
  
 * keep simStudy estimate stdErr parameter;
  
 run; proc sort; by simStudy parameter; run;


 proc datasets library=work noprint;
  delete __tempData__ __tempHist__;
 run;
 quit; 
 
%mend fit_Approx;