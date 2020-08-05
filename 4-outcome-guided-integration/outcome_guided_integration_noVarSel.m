%% Outcome-guided integration of the pancancer datasets of Hoadley et al.
close all
clc

cd '/Users/alessandracabassi/OneDrive - University of Cambridge, MRC Biostatistics Unit/PHD-PROJECTS/psms-pancancer-analysis/4-outcome-guided-integration'
addpath '/Users/alessandracabassi/OneDrive - University of Cambridge, MRC Biostatistics Unit/PHD-PROJECTS/SimpleMKL'
addpath '/Users/alessandracabassi/OneDrive - University of Cambridge, MRC Biostatistics Unit/PHD-PROJECTS/svm-km'

folder_path = '/Users/alessandracabassi/OneDrive - University of Cambridge, MRC Biostatistics Unit/PHD-PROJECTS/psms-pancancer-analysis/';
hpc_folder_path = '/home/ac2051/rds/hpc-work/psms-pancancer-analysis/';

%% Initialise variables

n_layers = 4; 
N = 2421;

%% Load PSMs

all_psms = zeros(N,N,n_layers);

% Copy number data

load(strcat(folder_path, 'mdi/psms-mat/psm_CN_average__noCuda.mat'))
all_psms(:,:,1) = average_psm;

% Methylation data
load(strcat(folder_path, 'mdi/psms-mat/psm_methylation_average__noCuda.mat'))
all_psms(:,:,2) = average_psm;

% miRNA data
load(strcat(folder_path, 'mdi/psms-mat/psm_miRNA_average__noCuda.mat'))
all_psms(:,:,3) = average_psm;

% Protein expression data
load(strcat(folder_path, 'mdi/psms-mat/psm_RPPA_average__noCuda.mat'))
all_psms(:,:,4) = average_psm;
%% Response

load(strcat(folder_path, 'data/tissue.mat'))
y = tissue;

%%  Initalize parameters of the algorithm 
% Parameters are similar to those used for mklsvm

C = 100;
lambda = 1e-7;
verbose = 1;

options.algo='oneagainstone';
options.seuildiffsigma=1e-4;
options.seuildiffconstraint=0.1;
options.seuildualitygap=1e-2;
options.goldensearch_deltmax=1e-1;
options.numericalprecision=1e-8;
options.stopvariation=1;
options.stopKKT=0;
options.stopdualitygap=1;
options.firstbasevariable='first';
options.nbitermax=500;
options.seuil=0.;
options.seuilitermax=10;
options.lambdareg = 1e-6;
options.miniter=0;
options.verbosesvm=0;
options.efficientkernel=0;

%%  Training

n_class = 10;
[beta,w,w0,pos,nbsv,~,~] = mklmulticlass(all_psms,y,C,n_class,options,verbose);
       
%% Write weights to file

save(strcat(folder_path, 'outcome-guided-integration-output/outcome_guided_integration_noVarSel_weights.mat'), 'beta');
