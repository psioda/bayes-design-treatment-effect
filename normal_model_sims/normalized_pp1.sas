***********************************************************************;
* Project           : Bayesian Clinical Trial Design Using Historical 
*                     Data that Inform the Treatment Effect
*
* Program name      : normalized_pp1.sas
*
* Author            : Matthew Psioda
*
* Date created      : 20171109
*
* Purpose           : Perform simulations for normalized power prior.
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
libname out "&root.\results\npp";

ods noresults;

** number of simulations per loop;
%let nSim   = 10000;

** number of loops;
%let nLoops = 20;

%put %upcase(no)TE: Total number of simulations = %eval(&nSim.*&nLoops.);

** identifier for historical dataset;
%let hd = 2;

** sample size in new trial;
%let N1 = 168;

** mean of beta prior for a0;
%let r0 = 0.41;

** dispersion parameter for beta prior for a0;
%let p0 = 5.00;

** simulation ID - will append to name of the results dataset;
%let sysparm = 1;

proc IML;
 ** seed;
call streaminit(15315);
call randseed(97954);

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
 N1 = &N1.;

 pi     = constant('pi');
 phi    = 0.975;


 r0 = &r0.;

 p0 = &p0.;

 ** points at which to approximate integrals with trapezoidal rule;
 ** this is used to integrate over a0 for posterior inference on gamma;
 ** the partition of [0,1] is very fine and therefore highly accurate;
 a0 = 0.0000000001||0.000000001||0.00000001||0.0000001||0.000001||0.00001||0.0001||
      do(0.0025,0.9975,0.0025)||
      0.999||0.9999||0.99999||0.999999||0.9999999||0.99999999||0.999999999||0.9999999999;

 ** compute a0 prior on log scale;
 logPrior_a0 = logPDF('beta',a0,r0*p0,(1-r0)*p0);

 ** construct intervals for a0 for trapezoidal rule;
    lid = 1:(ncol(a0)-1);
    hid = 2:ncol(a0);

    low   = a0[ lid ];
    high  = a0[ hid ];


  ** read the historical data; 
  use hist.hist_&hd.;
    read all var {y} into y0;
  close hist.hist_&hd.;

  ** compute summary statistics;
  y0Bar = y0[:];
  y0SS  = (y0##2)[+];
  SSE0  = ((y0 - y0Bar)##2)[+];
  MSE0  = SSE0/N0;

  pp = cdf('normal',0,y0Bar,sqrt(vr/N0));


  ** find support bounds for the truncated sampling priors;;
   den1 = pdf('normal',y0Bar,y0Bar,sqrt(vr/N0));
   ratio=1;
   xValue = y0Bar;
   do until(ratio>2);
    xValue = xValue + 1e-6;
     den2 = pdf('normal',xValue,y0Bar,sqrt(vr/N0));
     ratio = den1/den2;
   end;

   delta = xValue - y0Bar;

   lowEff  = min(y0Bar+delta,0);
   highEff = y0Bar-delta;


   call symput('lowEff',char( lowEff  ));
   call symput('HighEff',char( highEff ));


   den3 = pdf('normal',0,y0Bar,sqrt(vr/N0));
   ratio=1;
   xValue = 0;
   do until(ratio>2);
    xValue = xValue + 1e-6;
     den4 = pdf('normal',xValue,y0Bar,sqrt(vr/N0));
     ratio = den3/den4;
   end;

   lowNull = xValue;

   call symput('lowNull',char(lowNull));



 do replication = 1 to nLoops;


     ** obtain null and alternative tuncated sampling priors;
    NSP = randTN(nSim, y0Bar, sqrt(vr/N0),  0.00, &lowNull.);
    ASP = randTN(nSim, y0Bar, sqrt(vr/N0), &HighEff., &lowEff.);

 do scenario = 1 to 4;
   if scenario = 1 then SP = NSP;
   if scenario = 2 then SP = NSP * 0;
   if scenario = 3 then SP = ASP;
   if scenario = 4 then SP = ASP * 0 + y0Bar;

 do j = 1 to nSim;

   ** simulate newtrial data;
   yDat   = rand('normal',J(N1,1,SP[j]),sqrt(vr));
   yMean  = yDat[:];
   ySS    = (yDat##2)[+];
   ySSE   = (yDat##2)[+] - N1*yMean**2;
   yMSE   = ySSE/N1;
   vrD    = vr/N1;


   ** prior mean and variance;
   vr0     = vr/(N0*r0);
   mu0     = y0Bar;

   ** posterior mean and variance;
   vr1     = 1/(1/vrD+1/vr0);
   mu1     = (yMean/vrD+mu0/vr0)*vr1;

   ** posterior probability;
   postPostFixed  = 1 - cdf('normal',mu1/sqrt(vr1),0,1);


   ** compute posterior mean conditional on a0;
   vr0 = vr/(N0*a0);
   mu0 = y0Bar;
 
   vr1 = 1 / (1/vr0 + 1/vrD);
   mu1 = (yMean/vrD + mu0/vr0)#vr1;

   tSS  = ySS + a0*y0SS;
   tN   = N1 + N0*a0; 
   tMSE = (tSS - tN#mu1##2) / tN;

   ** compute posterior density for a0 evaluated at each a0;
   logPost_a0 = logPrior_a0 +
                - 0.5*log(2*pi*vr0) + 0.5/vr0*MSE0 + 0.5*(N0*a0)*log(2*pi*vr) + /** NPP normalizing constant **/
                + 0.5*log(2*pi*vr1) - 0.5/vr1#tMSE - 0.5*(tN   )*log(2*pi*vr) ;
   logPost_a0 = logPost_a0 - max(logPost_a0);
   post_a0 = exp(logPost_a0);

  ** compute normalizing constant using trapezoidal rule;
  low_p   = post_a0[ lid ];
  high_p  = post_a0[ hid ];
  area    = (0.5*(high-low)#(high_p + low_p))[+];

  ** normalize the posterior;
  Post_a0 = Post_a0 / area;

   ** compute prior probability of treatment effect * density weight
      over grid of a0 values;;
   PostProbWeight = cdf('normal',0,mu1,sqrt(vr1))#Post_a0;

   low_p  = PostProbWeight[ lid ];
   high_p = PostProbWeight[ hid ];

   ** compute posterior probability with trapezoidal rule;
   PostProb = (0.5*(high-low)#(high_p + low_p))[+];

   ** compute posterior mean for a0;
   ** compute prior probability of treatment effect * density weight;
   PostMeanWeight = a0#Post_a0;
   low_p  = PostMeanWeight[ lid ];
   high_p = PostMeanWeight[ hid ];

   ** compute posterior probability with trapezoidal rule;
   PostMean = (0.5*(high-low)#(high_p + low_p))[+];
 

   if j = 1 & scenario = 1 then do;
    tempYDat = yDat[,1]||repeat(PostProb,nrow(yDat),1)||repeat(postMean,nrow(yDat),1);
     create simDat from tempYDat[c={"y" "PostProb_Efficacy" "PostMean_a0"}];
      append from tempYDat;
     close simDat;
   end;

   if j > 1 then do;
     postProbAll = postProbAll // postProb;
      PostMeanAll = PostMeanAll // PostMean;
      postPostFixedAll = postPostFixedAll // postPostFixed;
   end;
   else do;
      postProbAll = postProb;
       PostMeanAll = PostMean;
       postPostFixedAll = postPostFixed;
   end;

 end;

  rejRate  = (postProbAll>phi)[:];
  rejRateFixed = (postPostFixedAll>phi)[:];

  PostMean = median(PostMeanAll);
  results = results // (&hd.||replication||scenario||r0||p0||rejRate||PostMean||rejRateFixed);

 end;
 end;



 create results from results[c={ "hist" "replication" "scenario" "r0" "p0" "rejRate" "PostMean_a0" "rejRateFixed_a0" }];
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

proc means data = results nway noprint;
 class hist scenario SamplingPrior r0 p0 / missing;
 var rejRate PostMean_a0 rejRateFixed_a0 ;
  output out = results mean= / autoname;
run;

data out.npp_%sysfunc(putn(&sysparm,z4.));
 set results;
 drop _type_ _freq_;
run;




ods rtf file = "&root.\results\npp\checkIMLProgramming.rtf";

proc univariate data = hist.hist_&hd. noprint;
 var y;
 output out = hist_stat n=n0 mean=y0bar var=s2;
run;

data _null_;
 set hist_stat;
 call symput('s2',strip(put(s2,best.)));
 call symput('n0',strip(put(n0,best.)));
run;

data total;
 set hist.hist_&hd.(in=a) simDat(keep=y);
  if a then hist = 1;
  else hist = 0;
   if _n_ = 1 then set simDat(drop=y);
run;


** mcmc code to check example;
 proc mcmc data = total nmc=50000 monitor=(gamma sig a0) plots=(none);
  by PostProb_Efficacy PostMean_a0;

  parm gamma 0 / slice;
  parm a0 0.5 / slice; 


  beginnodata;
  logPrior = logPDF('beta',a0,&r0.*&p0.,(1-&r0.)*&p0.)
           + &n0.*a0*0.5*log(2*constant('pi')) - 0.5*log(2*constant('pi') / (&n0.*a0))
           + 0.5*a0*(&n0.-1)*&s2.;
  endnodata;

  prior a0 ~ general(logPrior);
  prior gamma ~ general(0);

  if hist = 1 then weight = a0;
  else weight = 1.0;

  logModel = weight*logPDF('normal',y,gamma,1);

  sig = (gamma<0);

  model y ~ general(logModel);
run;
quit;

ods rtf close;
