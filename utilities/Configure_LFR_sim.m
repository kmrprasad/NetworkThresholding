thresholdpath = [pwd '/threshold_analysis'];
addpath(thresholdpath);
NMIpath=[pwd '/NMI'];
addpath(NMIpath);
BCTpath = [pwd '/BCT'];
addpath(BCTpath);
DATApath = [pwd '/data'];
addpath(DATApath);

data=[DATApath '/LFR_sim_data_20220622.mat']
load(data);

clear utilitiesPath NMIpath BCTpath data DATApath
