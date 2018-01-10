# bayes-design-treatment-effect
SAS programs used for Biostatistics manuscript "Bayesian Clinical Trial Design Using Historical Data that Inform the Treatment Effect".


==============
PART A: Design Application with Proportional Hazards Model 
==============

====
SAS Programs and Macros 
====
The following programs can be used to reproduce the simulation results for the case study application in the paper. The above programs are named 
in the order they MUST be run. The output from one step is used as input for the next step. Each of the SAS programs has a header that describes the
purpose of the program in the overall design process. If the programs are run as given (apart from updating paths to folders), they should reproduce the
results presented in the manuscript exactly. Version 9.4 of the SAS Software was used when writting these programs.

[0] 00_process_e1684_data.sas        --> Process data for e1684 study to create subject-level SAS dataset;

[1] 01_simulation_controls_null.sas  --> Create a dataset of simulation settings and random number seeds to be used in data generation and design evaluation for stage one of the simulation process
                                         where the focus is on computing the Bayesian type I error rate for a large number of combinations of target numbers of events and values of a0; The following
										 print out of the first few observations in the simulation settings dataset produced by this program should give users the main idea behind this dataset.
										 
																  Simulations        Sampling         Nuisance      Target              Data       Parameter
										  Compute   Simulation       Per              Prior          Parameters     Number    Value   Simulation    Sampling
											Node     Iteration       Loop            Dataset           Shared     of Events   of a0      Seed         Seed

											 1           1           500       sPrior.default_null        0          280        0       84150228    167955292
											 1           2           500       sPrior.default_null        0          280        0       31086397   1046511770
											 1           3           500       sPrior.default_null        0          280        0      737619331    980428618
											 1           4           500       sPrior.default_null        0          280        0      545235228    984462409
											 1           5           500       sPrior.default_null        0          280        0      706712105    394250946
											 1           6           500       sPrior.default_null        0          280        0      845876682    909962645
											 1           7           500       sPrior.default_null        0          280        0      804189096    461047745
											 1           8           500       sPrior.default_null        0          280        0      340793576    190291099
											 1           9           500       sPrior.default_null        0          280        0       71059337    265082720
											 1          10           500       sPrior.default_null        0          280        0      765919323    479432657
										 
										 Each observation in the dataset defines a set of possible design parameters. The first row in this example above indicates that the first compute node
										 is going to simulate 500 hypothetical new trials using the default null sampling prior. Information will only be borrowed through the treatment effect.
										 The new trial will target 280 events and the value of a0 will be 0. The last two numbers are seeds needed for reproducibility. In SAS (unfortunately), the
										 random number stream must be initialized again in each DATA step and/or PROC step.
										 
[2] 02_sampling_priors.sas           --> Reduce subject-level historical dataset to a dataset of sufficient statistics and create discrete approximations for all sampling priors;

[3] 03_calibrate_a0.sas              --> (To be run multiple times in parallel, if possible) Estimate the Bayesian type I error rate for a given set of design parameters (see [1] above)
										 
[4] 04_simulation_controls_alt.sas   --> Process the results from the first stage of the simulation to identify the optimal value of a0 for each target number of events and null sampling prior considered.
                                         Create a dataset of simulation settings and random number seeds to be used in data generation and design evaluation for stage two of the simulation process
                                         where the focus is on computing Bayesian power for valid combinations of the target number of events and values of a0. The following
										 print out of the first few observations in the simulation settings dataset produced by this program should give users the main idea behind this dataset.

																				  Null            Alternative
															   Simulations      Sampling            Sampling       Nuisance    Target             Data     Parameter
											Compute Simulation     Per            Prior              Prior        Parameters   Number   Value  Simulation  Sampling
											  Node   Iteration     Loop          Dataset            Dataset         Shared   of Events  of a0     Seed       Seed

											   1         1         500     sPrior.default_null sPrior.default_alt      0        280    0.90727   78904831  969015141
											   1         2         500     sPrior.default_null sPrior.default_alt      0        280    0.90727 1045488691   54992238
											   1         3         500     sPrior.default_null sPrior.default_alt      0        280    0.90727 1034750852  240510388
											   1         4         500     sPrior.default_null sPrior.default_alt      0        280    0.90727  268996136   29726410
											   1         5         500     sPrior.default_null sPrior.default_alt      0        280    0.90727  531349376  819324977
											   1         6         500     sPrior.default_null sPrior.default_alt      0        280    0.90727  935254060  381897243
											   1         7         500     sPrior.default_null sPrior.default_alt      0        280    0.90727  568822962  389770595
											   1         8         500     sPrior.default_null sPrior.default_alt      0        280    0.90727  714879593 1004695207
											   1         9         500     sPrior.default_null sPrior.default_alt      0        280    0.90727  353048913 1025541179
											   1        10         500     sPrior.default_null sPrior.default_alt      0        280    0.90727  915764804  712061711

										 This dataset is similar in structure and purpose to that discussed in [1].
										 
[5] 05_calculate_power.sas           --> (To be run multiple times in parallel fashion, if possible) Estimate the Bayesian power for a given set of design parameters (see [4] above).
										 
[6] 06_process_power.sas             --> Process the results from the second stage of the simulation process to identify the smallest number of events that yields appropriate Bayesian power given the
                                         Bayesian type I error restriction. Create a dataset of the optimal design parameters for each null sampling prior and alternative sampling prior combination.
										 After this step is complete users will have a SAS dataset containing parameters of the optimal designs. The following print out describes the contents.
										 
												   Null                 Alternative
												 Sampling                Sampling            Nuisance       Target                Estimated
												   Prior                   Prior            Parameters      Number      Value      Bayesian
												  Dataset                 Dataset             Shared      of Events     of a0       Power

											sPrior.default_null    sPrior.default_alt            0           525       1.00000     0.79975
											sPrior.default_null    sPrior.freq_null              0           355       0.99686     0.07084
											sPrior.default_null    sPrior.general_alt            0           515       1.00000     0.80029
											sPrior.default_null    sPrior.point_mass_alt         0           305       0.94661     0.79874
											sPrior.freq_null       sPrior.default_alt            0           815       0.00716     0.79951
											sPrior.freq_null       sPrior.freq_null              0           855       0.00868     0.02538
											sPrior.freq_null       sPrior.general_alt            0           785       0.00601     0.79981
											sPrior.freq_null       sPrior.point_mass_alt         0           445       0.00303     0.79881
											sPrior.general_null    sPrior.default_alt            0           575       0.83673     0.79983
											sPrior.general_null    sPrior.freq_null              0           745       0.99160     0.05870
											sPrior.general_null    sPrior.general_alt            0           565       0.82537     0.80006
											sPrior.general_null    sPrior.point_mass_alt         0           355       0.60065     0.79567
										 
										 
										 

The SAS programs above make use of the following SAS macros which are described briefly:

[1] collapse.sas       --> This SAS macro is used to collapse a subject-level time-to-event dataset down to the sufficient statistics for the 
                           stratified proportional hazards model. The sufficient statistics are the number of events and total time at risk 
                           associated with each value of the hazard ratio covariates (including treatment indicator) and baseline hazard component. 
						   
[2] fix_approx.sas     --> This SAS macro is used fit the PH model using the power prior via weighted maximum likelihood approximation using PROC GENMOD. 

[3] fit_mcmc.sas       --> This SAS macro is used fit the PH model using the power prior via MCMC using PROC MCMC.

[4] simulateStudy.sas  --> This SAS macro is used to simulate data for the new trial. The macro simulates data from a proportional hazards model with piecewise constant
                           baseline hazard.

The SAS macro programs contain header documentation that gives more detail on the purpose of the macro and the macro variables required by the macro. 


====
Directory Structure
====

The following directory structure is assumed for the tools (as written) and a brief description of the contents of folders is given:

[root]/data

[root]/data/raw_data                             --> contains raw data (text file) for the e1684 study that is used as input for the design (i.e., the historical data). The SAS programs above create several SAS datasets
                                                     and store them in this folder. Only the raw text file is included in the GitHub repository to keep the size down. All other content can be generated by the SAS programs above.

[root]/data/sampling_priors                      --> contains SAS datasets that correspond to the discrete approximations to the various sampling priors used in the paper. (empty on GitHub)
 
[root]/data/simulation_controls                  --> contains SAS datasets that correspond to the simulations settings used for stage one and stage two simulations. These datasets are used as input
                                                     for programs [3] and [5] above. (empty on GitHub)

[root]/data/simulation_results             
      
[root]/data/simulation_results/calibrate_a0      --> contains the analysis results (e.g. estimated Bayesian type I error rate) for designs corresponding to the simulation settings investigated in stage one. (empty on GitHub)

[root]/data/simulation_results/calculate_power   --> contains the analysis results (e.g. estimated Bayesian power) for designs corresponding to the simulation settings investigated in stage two. (empty on GitHub)

[root]/macros                                    --> contains the SAS macros described above.

[root]/cluster-scripts                           --> Shell scripts to run jobs on a Linux cluster using a SLURM scheduler. There is one batch script for each program. Programs [3] and [5]
                                                     are designed to be run as array jobs.

[root]/cluster-logs                              --> SAS logs for the program runs corresponding to the simulations presented in the paper.

[root]/cluster-out                               --> OUT files for produced by SLURM scheduler.

[root]/cluster-err                               --> ERR files for produced by SLURM scheduler.

[root]/                                          --> location where the SAS programs above reside.

In each program there is a macro variable named ROOT that stores the users root directory path. That path must be changed to match the users chosen directory for the SAS programs to work.

======
Parallel Execution of Programs
======

Programs [3] and [5] above are designed to be run as array jobs on a computing cluster to perform simulation studies over a grid of design inputs in parallel fashion. For
example, the authors use the following array script for a SLURM scheduler on a computing cluster for program [3]:

############ Begin Script #########

#!/bin/bash

#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 48:00:00
#SBATCH --mem 4000
#SBATCH --output=./../cluster-out/03-%a.out
#SBATCH --error=./../cluster-err/03-%a.err
#SBATCH --array=1-600

## add SAS
module add sas/9.4


## run SAS command
sas -work /dev/shm -noterminal ./../03_calibrate_a0.sas -log "./../cluster-logs/03-$SLURM_ARRAY_TASK_ID.log" -sysparm "$SLURM_ARRAY_TASK_ID"

############ End Script #########

Users can always run programs [3] and [5] repeatedly on a single computer while cycling through values for SYSPARM until all design inputs have been explored. For each run of the program,
the datasets created will use the value of SYSPARM as a suffix for the dataset names (e.g., calibrate_a0_00001.sas7bdat for SYSPARM=1). Note that the number of times
programs [3] and [5] need to be run (each using a different value of SYSPARM) will depend on the number of design inputs to be investigated. For the paper, we performed a very large
number of simulations (different null sampling priors, alternative sampling priors, many possible number of events, many possible a0 values) and so using a cluster to run the programs many
times in parallel was useful.

==============
PART B: Design Simulations with Simple Normal Model 
==============

====
Directory Structure
====

All programs for the normal model simulation expect the following directory structure: 

[root]/hist-dat       - Contains the 5 hypothetical historical datasets of varying levels of informativeness
[root]/results
[root]/results/amap   - Contains results of optimal adaptive design using robust mixture prior (sample size adaptively increased)
[root]/results/map    - Contains results of optimal design using robust mixture prior 
[root]/results/npp    - Contains results of optimal design using the normalized power prior
[root]/results/pp     - Contains results of optimal design using power prior, supervised power prior, and modified critical value approach

In each program there is a macro variable named ROOT that stores the users root directory path. That path must be changed to match the users chosen directory for the SAS programs to work.

====
Historical Datasets
====
Simulated Historical Datasets
All 5 stimulated historical datasets used in the paper are included in the subfolder named “HIST-DAT”. Users should note that the historical datasets are 
named “HIST_X” where X is an integer from 1 to 5. Please be aware of the following mapping between datasets and informativeness of the data;

HIST_1 - posterior probability that $\gamma<0$ = 0.975
HIST_2 - posterior probability that $\gamma<0$ = 0.950
HIST_3 - posterior probability that $\gamma<0$ = 0.900
HIST_4 - posterior probability that $\gamma<0$ = 0.850
HIST_5 - posterior probability that $\gamma<0$ = 0.999

Thus, the datasets are not named in descending order of informativeness.


====
Description of SAS Programs
====
SAS Programs that reproduce the simulation results presented in the manuscript are provided in the repo and are listed below:

[1] PP1.SAS - Estimates the Bayesian type I error rate and power using truncated sampling priors (K=2) and based on the 
              partial-borrowing power prior, the supervised power prior, and the modified critical value approach. In the body 
			  of the program, the user can specify design parameters that are used in the simulation. In particular, the 
			  following %LET statements define key macro variables:
			  
** identifier for historical dataset;
%let hd = 2;

** sample size in new trial;
%let N1 = 168;

** fixed value of a0 for partial-borrowing power prior;
%let a0 = 0.44;

** calibration parameter for supervised power prior;
%let s0 = 0.96;

** modified critical value;
%let phi0 = 0.955;

The current settings correspond to an optimal design for the simulated historical dataset with posterior probability of treatment 
efficacy equal to 0.95. The program is thoroughly commented to explain the various components and random number seeds can be 
specified in the PROC IML step. All other programs below incorporate random number seeds and are similarly commented.

[2] NORMALIZED_PP1.SAS - Estimates the Bayesian type I error rate and power using truncated sampling priors (K=2) and 
                         based on the normalized power prior. In the body of the program, the user can specify design parameters 
						 that are used in the simulation. In particular, the following %LET statements define key macro variables:
						 
** identifier for historical dataset;
%let hd = 2;

** sample size in new trial;
%let N1 = 168;

** mean of beta prior for a0;
%let r0 = 0.41;

** dispersion parameter for beta prior for a0;
%let p0 = 5.00;

This program takes significantly longer to run than others due to the computational burden of integrating out $a_0$ in order to 
compute the posterior probability that $\gamma<0$. Still, it only takes a few minutes to complete a highly accurate set of simulations. The 
final section of the program verifies the correctness of our IML implementation of the normalized power prior by comparing to results from 
the SAS MCMC procedure.

[3] MAP_PRIOR1.SAS - Estimates the Bayesian type I error rate and power using truncated sampling priors (K=2) and based on the two-component 
                     robust mixture prior. In the body of the program, the user can specify design parameters that are used in the simulation. 
					 In particular, the following %LET statements define key macro variables:
					 
** identifier for historical dataset;
%let hd = 2;

** sample size in new trial;
%let N1 = 168;

** prior weight for historical data;
%let weight = 0.74;

** variance inflation factor;
%let varFact = 1.0;

The final section of the program verifies the correctness of our IML implementation of the robust two-component mixture prior by comparing to 
results from the SAS MCMC procedure.

[4] AMAP_PRIOR1.SAS - Estimates the Bayesian type I error rate and power using truncated sampling priors (K=2) and based on the two-component 
robust mixture prior in an adaptive design where the sample size in the new trial can be adaptively increased. In the body of the program, 
the user can specify design parameters that are used in the simulation. In particular, the following %LET statements define key macro variables:

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


