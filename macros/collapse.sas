/*
This SAS macro is used to collapse a subject-level time-to-event dataset down to the sufficient statistics for the 
stratified proportional hazards model. The sufficient statistics are the number of events and total time at risk 
associated with each value of the hazard ratio covariates (including treatment indicator) and baseline hazard 
component. 

This macro requires the following parameters:
 
  MACRO VARIABLES REQUIRED:
   [1] fitHazComp   = list of the number of baseline hazard components in each stratum - each stratum separated by a pipe;
   [2] fitHazBreaks = list of the change points for baseline hazard in each stratum - each stratum separated by a pipe;
   [3] stratlevels  = number of strata;
   [4] dsIN         = dataset to be used as input for data reduction;
   [5] dsOUT        = dataset to be written out after data reduction;
   [6] writeBreaks  = indicator of whether or not a dataset should be generated to store the baseline hazard change points; 
                      (this is only needed once for the historical dataset)
   [7] dsBreaks     = dataset to be created that stores the baseline hazard change points;
                      (this is only relevant when writeBreaks = 1 but still needs to be defined in all cases)
*/

%macro collapse;
%let maxComp = 0;
%do s = 1 %to &stratLevels.;
  %if &maxComp < %scan(&fitHazComp.,&s,|) %then %let maxcomp = %scan(&fitHazComp.,&s,|);
%end; 


** create hazard breaks dataset;
data __breaks__;

 array t(%eval(&maxComp.+1)) %do k = 1 %to %eval(&maxComp.+1); tpt%eval(&k-1) %end;;
 %do s = 1 %to &stratLevels.;
   %let nComp  = %scan(&fitHazComp.,&s.,|);
   %let tParms = %scan(&fitHazBreaks.,&s.,|);


   riskStratum = &s.;
   nStratComp   = &nComp.;


   t(1) = 0;
   %do k = 1 %to %eval(&nComp.-1);
    t(%eval(&k.+1)) = %scan(&tParms.,&k.,%str( ));
   %end;
   t(%eval(&k.+1)) = 1e10;
   output;
 %end;
run;


 data __tempData__;
  set &dsIn.;
 run;

 proc sql undo_policy=none;
  create table __tempData__ as
  select a.*,b.nStratComp %do k = 1 %to %eval(&maxComp.+1); ,b.tpt%eval(&k-1) %end;
  from __tempData__ as a, __breaks__ as b
  where a.riskStratum=b.riskStratum;*
  order by a.simstudy,a.riskStratum,a.x1,a.obsTime;
 quit;

 data __tempData__;
  set __tempData__;
   array tpt(*) tpt:;
   do i = 1 to nStratComp;

    r        = 0;
	v        = 0;
    interval = i;

    if tpt(i) < obsTime <= tpt(i+1) then do;
      r = obsTime - tpt(i);
      v = 1-censor;
      i = nStratComp + 1;
	  output;
	end;
	else if obsTime > tpt(i+1) then do;
      r = tpt(i+1) - tpt(i);
      v = 0;
	  output;
	end;
   end;
run;

proc means data = __tempData__ noprint nway;
  class simstudy riskStratum x1 interval;
  var r v;
  output out = &dsOut.(drop=_type_ _freq_) sum=riskTime numEvents;
run;

 %if &writeBreaks = 1 %then %do;
 data &dsBreaks.;

   %do s = 1 %to &StratLevels.;
   %do j = 1 %to &maxComp.; 
     t&s.&j.= . ;
   %end;
   %end;

    %do s = 1 %to &StratLevels.;
	%let nComp = %scan(&fitHazComp.,&s,|);
	%let bks   = %scan(&fitHazBreaks.,&s,|);
	%do j = 1 %to %eval(&nComp.-1); t&s.&j. = %scan(&bks.,&j., %str( )); %end;
	%end;
run;
%end;

 proc datasets library=work noprint;
  delete __breaks__ __tempData__;
 run;
 quit; 


%mend collapse;