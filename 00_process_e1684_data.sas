/*
Program Name: 00_process_e1684_data.sas
Program Creator: Matthew Psioda;

This SAS program reads in the raw data file (named mina-e1684-e1690.txt), subsets down to the E1684 that was 
used for the case study presented in the paper, and then creates a subject-level SAS dataset that is stored 
in the folder: &root.\data\raw_data (RAWOUT SAS library)

*/

%let root        = /nas/longleaf/home/psioda/stat-projects/bayesDesignTreatment;
%let rawDataPath = &root./data/raw_data;

libname rawout "&rawDataPath.";

data rawout.e1684;
 infile "&rawDataPath/mina-e1684-e1690.txt" dlm=' ' missover firstobs=2;
 input case study age trt sex perform nodes breslow stage failtime rfscens survtime scens;

  ** one subject had a failure time of zero and this 
     code sets the value to something negligible for processing;
  if failtime = 0 then failtime=1e-5;

  ** remove e1690 data;
  if study = 1684;
run;

  ** the example in the paper is based on a complete case analysis
     restricted to stage 4 melanoma patients;
data rawout.e1684_stage4;
 set rawout.e1684;
 where stage = 4 and nodes>.;

 if nodes <= 2 then riskStratum = 1;
 else if nodes > 2 then riskStratum = 2;

 simStudy = 1;
 censor = 1-rfscens;
 keep trt nodes failtime censor riskStratum simStudy;
  rename failtime = obsTime trt = x1 ;
run;
