***********************************************************************;
* Project           : Bayesian Clinical Trial Design Using Historical 
*                     Data that Inform the Treatment Effect
*
* Program name      : amap_prior1.sas
*
* Author            : Matthew Psioda
*
* Date created      : 20171109
*
* Purpose           : Simulations for robust MAP prior (Adaptive).
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
libname out "&root.\results\amap";

** number of simulations per loop;
%let nSim   = 50000;

** number of loops;
%let nLoops = 20;

%put %upcase(no)TE: Total number of simulations = %eval(&nSim.*&nLoops.);

** identifier for historical dataset;
%let hd = 2;

** sample size in stage 1 of new trial;
%let N1 = 160;

** sample size in stage 2 of new trial;
%let N0Adp = 20;

** simulation ID - will append to name of the results dataset;
%let sysparm = 1;

** prior weight for historical data;
%let weight = 0.69;

** variance inflation factor;
%let varFact = 1.0;

proc IML;
 ** seed;
call streaminit(4362);
call randseed(7959);

** number of simulations to estimate properties;
nSim   = &nSim.;
nLoops = &nLoops.;

** sample from truncated normal;
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
 phiMod = 0.958;
 
 ** Increase in sample size due to adapting;
 N0Adp = &N0Adp.;

 ** initial weights for aMAP prior components;
 omega  = &weight.;
 omega  = omega ||(1-omega);


  ** read the historical data;  
  use hist.hist_&hd.;
    read all var {y} into y0;
  close hist.hist_&hd.;

  ** compute summary statistics;
  y0Bar = y0[:];
  y0Bar = y0[:];
  SSE0  = ((y0 - y0Bar)##2)[+];
  pp = cdf('normal',0,y0Bar,sqrt(vr/N0));


  ** parameters to define MAP prior components;
  mu0Vec = y0Bar || 0.000;
  vrFrac = 1.000 || &varFact.;

  varFact = &varFact.;

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

   call symput('sp',strip(char(scenario)));

   ** simulate newtrial data;
   yDat   = rand('normal',repeat(SP,1,N1),sqrt(vr));

   yMean  = yDat[,:];
   ySSE   = (yDat##2)[,+] - N1*yMean##2;
   vrD    = vr/N1;
 

   ** analysis based on MAP prior (non-adaptive);
   weight = repeat(log(omega),nSim,1);
   pp     = J(nSim,2,0);

   do j = 1 to 2;

      ** prior mean and variance;
      vr0     = vr/N0*vrFrac[j]; 
      mu0     = mu0Vec[j];

       ** posterior mean and variance;
      vr1     = 1/(1/vrD+1/vr0);
      mu1     = (yMean/vrD+mu0/vr0)*vr1;

      weight[,j] = weight[,j] - 0.5*( ySSE/vrD + mu0**2/vr0 - mu1##2/vr1 )
                  - N1*0.5*log(2*pi*vr) - 0.5*log(2*pi*vr0) + 0.5*log(2*pi*vr1);

      pp[,j] = 1 - cdf('normal',mu1/sqrt(vr1),0,1);

   end;

   weight  = exp(weight - weight[,<>]);
   weight  = weight / weight[,+];

   ** Adaptive MAP Design;

     ** Determine when studies adapt;
     ** adapt when the posterior weight is less than starting weight;
     adapt   = (weight[,1] < omega[1]);

     ** simulate newtrial data for N0Adp more patients;;
     yAdpDat   = rand('normal',repeat(SP,1,N0Adp),sqrt(vr));

     yDatNew = yDat || (yAdpDat#(adapt=1));
     SSAdp   = (N1 + N0Adp*(adapt=1));

     yAdpMean  = yDatNew[,+] / SSAdp;
     yAdpSSE   = (yDatNew##2)[,+] - SSAdp#yAdpMean##2;
     vrAdp     = vr/SSAdp;

        ** for each value of simulation study, compute the new weights;
        weightAdp = repeat(log(omega),nSim,1);
        ppAdp     = J(nSim,2,0);
        do j = 1 to 2;
           ** prior mean and variance;
           vr0     = vr/N0*vrFrac[j]; 
           mu0     = mu0Vec[j];

           ** posterior mean and variance;
           vr1     = 1/(1/vrAdp+1/vr0);
           mu1     = (yAdpMean/vrAdp+mu0/vr0)#vr1;

           weightAdp[,j] = weightAdp[,j] - 0.5*( yAdpSSE/vrAdp + mu0**2/vr0 - mu1##2/vr1 )
                       - SSAdp#0.5*log(2*pi*vr) - 0.5*log(2*pi*vr0) + 0.5*log(2*pi*vr1);

           ppAdp[,j] = 1 - cdf('normal',mu1/sqrt(vr1),0,1);

        end;

        weightAdp  = exp(weightAdp - weightAdp[,<>]);
        weightAdp  = weightAdp / weightAdp[,+];
        rejRateAdp = ((weightAdp#ppAdp)[,+] > phi)[:];
        AdtPrb     = adapt[:];
          


  medAdp_weight = median(weightAdp[,1]);
  meanAdp_weight = mean(weightAdp[,1]);

  ssAdp = ssAdp[:];


  results = results // (&hd.||replication||scenario||varFact||N0Adp||N1||
                                 omega[1]||AdtPrb||medAdp_weight[1]||meanAdp_weight[1]||ssAdp||rejRateAdp);

 end;
 end;

 create results from results[c={ "hist" "replication" "scenario" "varFact" "N0Adp" "N1"
                                "priorwAdp" "prbAdapt" "medwAdp" "meanwAdp" "ssAdp" "rejRateAdp" }];
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

proc means data = results noprint nway;
 class hist scenario SamplingPrior varFact N0Adp N1 priorwAdp / missing;
 var  prbAdapt medwAdp meanwAdp ssAdp rejRateAdp;
 output out = results mean =  prbAdapt medwAdp meanwAdp ssAdp rejRateAdp;
run;

data out.amap_%sysfunc(putn(&sysparm,z4.));
 set results;
 drop _type_ _freq_;
run;


