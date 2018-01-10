/*
Program Name: 01_simulation_controls_null.sas
Program Creator: Matthew Psioda;

This SAS program generates the first stage simulation settings that were explored in the SDD example based on the e1684 dataset that is 
presented in the paper. Initially, we needed to determine the maximum value of a0 for each number of events in the set under 
consideration (250 events to 850 events with a 10 event step size).

This program generates a dataset that contains an large array combinations of target number of events and a0 values in a SAS dataset and 
assigns them to different compute node IDs for parallel simulations on a cluster; The NODE_IDX variable identifies the compute 
node by number. Thus, if a user is going to use 600 compute nodes for computation as we did, NODE_IDX will take values between 1 
and 600. Using this number of compute nodes is not necessary but will result in fast completion of simulations if such 
resources are available.

The created SAS dataset of simulation settings is written to the location: &root./data/simulation_controls;
The dataset is named controls.simulation_controls_null and contains the following variables:
[1] node_idx – Integer value that identifies the compute node to be used for a given set of simulations;
[2] nPerLoop – number of datasets to simulate in a given loop cycle – this number cannot be too large due to system 
    resources. A value of 500 is used in the current setup. The more datasets that are simulated at once, the more 
	memory or disk space needs to be available to SAS to process them;
[3] sp – name of the dataset that stores the sampling prior results for the given simulation setting;
[4] shareBLH – indicator variable for whether or not all parameters are shared between model (must be zero in current software);
[5] nuTarget – desired number of events in the simulated new trials;
[6] a0 – value of a0 for the power prior;
[7] inner_idx – observation counter for simulation settings to be performed by a given compute node;
[8] seed – random number seed for new trial data generation;
[9] parmSampleSeed – random number seed for drawing samples from discrete approximation to the sampling priors;

The dataset created by this program is stored in the folder: &root./data/simulation_controls (CONTROLS SAS library)
*/

%let root     = /nas/longleaf/home/psioda/stat-projects/bayesDesignTreatment;
%let contPath = &root./data/simulation_controls;

proc datasets library=work noprint kill; run; quit;
option compress=Yes;


libname _all_ clear;
libname controls "&contPath.";


%let number_of_nodes = 600; ** number of compute nodes on cluster;
%let nPerLoop        = 500; ** number of simulated datasets per simulation loop;
%let nLoops          = 200; ** number of simulation loops;

** the design operating characteristics will be estimated based on &nPerLoop * &nLoops simulation studies;

data controls.simulation_controls_null;
 nPerLoop = &nPerLoop.;
 nLoops   = &nLoops.;
 idx      = 1;


 length sp $30.;
 do sp        = "sPrior.default_null","sPrior.freq_null","sPrior.general_null"; ** null sampling priors to explore;
 do shareBLH  = 0 to 0;                                                         ** indicator for whether or not nuisance parameters are shared;
 do nuTarget  = 280 to 850 by 10;                                               ** possible sample size (number of events);
 do a0        = 0.00,0.01,0.02,0.03,0.04,0.05,0.075 to 1.00 by 0.025;           ** possible values for a_0;


  output;  idx + 1; 

 end;
 end;
 end;
 end;

 call symput('num_idx',strip(put(idx-1,best.)));
 drop idx;
run; 
proc sort; by sp shareBLH nuTarget a0; run;
%put NO%upcase(te):There will be &=num_idx different scenarios considered in this simulation;

** this DATA step allocations the different simulations to the compute nodes;
data node_id;
 stop = 0;
 i=0;
 do outer_idx = 1 to 20000;
 do node_idx = 1 to &number_of_nodes.;
   i+1;
   output;
   if i = &num_idx. then stop;
 end;
 end;
 keep node_idx;
run; proc sort; by node_idx; run;

data controls.simulation_controls_null;
 merge node_id controls.simulation_controls_null ;
run;


data controls.simulation_controls_null(index=(node_idx));
 set controls.simulation_controls_null;
 call streaminit(123512);

 ** two random number seets are needed for each simulation in order
    to ensure reproducibility and they are generated here;
 do inner_idx = 1 to nLoops;
   seed           = round(1 + rand('uniform')*2**30);
   parmSampleSeed = round(1 + rand('uniform')*2**30);
    output;
 end;
 drop nLoops;
run;


/*  proc contents data = controls.simulation_controls_null ; run;  */
