***********************************************************************;
* Project           : Bayesian Clinical Trial Design Using Historical 
*                     Data that Inform the Treatment Effect
*
* Program name      : pp1.sas
*
* Author            : Matthew Psioda
*
* Date created      : 20171109
*
* Purpose           : Simulations for partial-borrowing power prior,
*                     supervised power prior, and modified critical value
*                     approach.
*
* Revision History  :
*
* Date        Author      Ref    Revision (Date in YYYYMMDD format) 
*
***********************************************************************;
** specify the root folder for the simulation files;
** the root folder needs subfolders named as follows:
   [1] .\hist-data
   [2] .\results
   [3] .\results\pp
   [4] .\results\map
   [5] .\results\amap
   [6] .\results\npp;

%let root = C:\Users\psioda\Documents\Research Papers Temp\bayesDesignTreatment\sas_programs_github\normal_model_sims;

ods _all_ close;

proc datasets library=work noprint kill; run; quit;

** library for simulated historical datasets;
libname hist "&root.\hist-dat";

** library for simulation results;
libname out "&root.\results\pp";

** number of simulations per loop;
%let nSim   = 50000;

** number of loops;
%let nLoops = 20;

%put %upcase(no)TE: Total number of simulations = %eval(&nSim.*&nLoops.);

** identifier for historical dataset;
%let hd = 2;

** sample size in new trial;
%let N1 = 168;

** simulation ID - will append to name of the results dataset;
%let sysparm = 1;

** fixed value of a0 for partial-borrowing power prior;
%let a0 = 0.44;

** calibration parameter for supervised power prior;
%let s0 = 0.96;

** modified critical value;
%let phi0 = 0.955;




proc IML;

 ** seed;
call streaminit(5737537);
call randseed(353515);

** number of simulations to estimate properties;
nSim   = &nSim.;
nLoops = &nLoops.;

** sample from truncated normal - code by Rick Wicklin;
** https://blogs.sas.com/content/iml/2013/07/24/the-truncated-normal-in-sas.html;
start randTN(N, mu, sigma, a, b);
   /* Support one-sided truncation. Define Phi(.M)=0 and Phi(.P)=1 */
   cdf_a = choose(a=., 0, cdf("Normal", a, mu, sigma));  /* Phi(a) */
   cdf_b = choose(b=., 1, cdf("Normal", b, mu, sigma));  /* Phi(b) */
   u = j(N,1,.);
   call randgen(u, "Uniform");           /* U ~ U(0,1) */
   z = cdf_a + (cdf_b - cdf_a)*u;        /* Z ~ U(Phi(a), Phi(b)) */
   x = quantile("Normal", z, mu, sigma); /* x always in [a,b] */
   return(x);
finish;


 ** known variance;
 vr = 1.00;

 ** historical sample size;
 N0 = 50;

 ** sample size in new and historical studies;
 N1   = &N1.;

 pi     = constant('pi');
 phi    = 0.975;

 ** modified critical value;
 phiMod = &phi0.;
 
 ** modified prior-data confict statistic;
 s0 = &s0.;

 ** borrowing parameter for power prior;
 a0  = &a0.;


  ** read the historical data; 
  use hist.hist_&hd.;
    read all var {y} into y0;
  close hist.hist_&hd.;

  ** compute summary statistics;
  y0Bar = y0[:];
  y0Bar = y0[:];
  SSE0  = ((y0 - y0Bar)##2)[+];
  pp = cdf('normal',0,y0Bar,sqrt(vr/N0));

  


  ** parameters needed to find support bounds 
     for the truncated sampling priors;;
  mu0Vec = y0Bar || 0.000;
  vrFrac = 1.000 || 1.00;


   den1 = pdf('normal',mu0Vec[1],mu0Vec[1],sqrt(vr/N0));
   ratio=1;
   xValue = mu0Vec[1];
   do until(ratio>2);
    xValue = xValue + 1e-6;
     den2 = pdf('normal',xValue,mu0Vec[1],sqrt(vr/N0));
     ratio = den1/den2;
   end;

   delta = xValue - mu0Vec[1];

   lowEff  = min(mu0Vec[1]+delta,0);
   highEff = mu0Vec[1]-delta;


   call symput('lowEff',char( lowEff  ));
   call symput('HighEff',char( highEff ));


   den3 = pdf('normal',0,mu0Vec[1],sqrt(vr/N0));
   ratio=1;
   xValue = 0;
   do until(ratio>2);
    xValue = xValue + 1e-6;
     den4 = pdf('normal',xValue,mu0Vec[1],sqrt(vr/N0));
     ratio = den3/den4;
   end;

   lowNull = xValue;
   call symput('lowNull',char(lowNull));

 
 do replication = 1 to nLoops;

     ** obtain null and alternative tuncated sampling priors;
    NSP = randTN(nSim, mu0Vec[1], sqrt(vr/N0),  0.00, &lowNull.);
    ASP = randTN(nSim, mu0Vec[1], sqrt(vr/N0), &HighEff., &lowEff.);


 do scenario = 1 to 4;

   if scenario = 1 then SP = NSP;
   if scenario = 2 then SP = NSP * 0;
   if scenario = 3 then SP = ASP;
   if scenario = 4 then SP = ASP * 0 + mu0Vec[1];

  ** call symput('sp',strip(char(scenario)));

   ** simulate newtrial data;
   yDat   = rand('normal',repeat(SP,1,N1),sqrt(vr));

   yMean  = yDat[,:];
   ySSE   = (yDat##2)[,+] - N1*yMean##2;
   vrD    = vr/N1;
 

** power prior analysis;

   ** prior mean and variance;
   vr0     = vr/(N0*a0);
   mu0     = mu0Vec[1];

   ** posterior mean and variance;
   vr1     = 1/(1/vrD+1/vr0);
   mu1     = (yMean/vrD+mu0/vr0)*vr1;

   ** posterior probability;
   pp  = 1 - cdf('normal',mu1/sqrt(vr1),0,1);
   rejRatePP = ((pp > phi) )[:];

 ** modified critical value analysis;

   ** posterior probability;
   ppMC      = 1 - cdf('normal',yMean/sqrt(vrD),0,1);
   rejRateMC = (ppMC > phiMod)[:];

 ** statistic-based borrowing;


   delta0  = exp(s0*(logpdf('normal',mu0,yMean,sqrt(vrD)) - logpdf('normal',yMean,yMean,sqrt(vrD))));

   ** prior mean and variance;
   vr0     = vr/(N0*delta0);
   mu0     = mu0Vec[1];

   ** posterior mean and variance;
   vr1     = 1/(1/vrD+1/vr0);
   mu1     = (yMean/vrD+mu0/vr0)#vr1;

   ** posterior probability;
   ppw  = 1 - cdf('normal',mu1/sqrt(vr1),0,1);
   rejRateStat = (ppw > phi)[:];

     delta0 = median(delta0);

  results = results // (&hd.||N1||replication||scenario||a0||rejRatePP||phiMod||rejRateMC||s0||delta0||rejRateStat);

 end;
 end;

 create results from results[c={ "hist" "ss" "replication" "scenario" "a0" "rejRatepp" 
                                "phiMod" "rejRateMC" "s0" "delta0" "rejRateStat" }];
  append from results;
 close results;

quit;

  proc format; 
   value fmtShort
    1 = 'TN'
    2 = 'PMN'
    3 = 'TA'
    4 = 'PMA';
  run;
 quit;

data results;
 set results;
  SamplingPrior = put(scenario,fmtShort.);
run;

** collapse results over loops;
proc means data = results noprint nway;
 class hist ss scenario SamplingPrior a0 s0 phiMod / missing;
 var rejRatepp rejRateMC delta0 rejRateStat;
 output out = results mean = rejRatepp rejRateMC delta0 rejRateStat;
run;



data out.pp_%sysfunc(putn(&sysparm,z4.));
 set results;
 drop _type_ _freq_;
run;



  

