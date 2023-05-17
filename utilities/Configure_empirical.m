% CHANGE TO YOUR SYSTEM FILE LOCATION
root = '/data1/NetworkThresholding/NetworkThresholding';

thresholding = sprintf('%s/thresholding',root);
addpath(thresholding);

utilities = sprintf('%s/utilities',root);
addpath(utilities);

BCT = sprintf('%s/BCT',root);
addpath(BCT);

NMI = sprintf('%s/NMI',root);
addpath(NMI);

visualization = sprintf('%s/visualization',root);
addpath(visualization);

data = sprintf('%s/empirical_UKB_HCP_100.mat',root);
load(data);

LFR = sprintf('%s/PowerLawNetMat.mat',root);
load(LFR);
modules_groundTruth = m';

clear thresholding utilities BCT NMI visualization  LFR 