%% load and configure
%no negatives


%configure
Configure_LFR_sim;

load('Results_noNegatives.mat');

CC_FC_noNeg = CC_FCsim_abs; %this one has a bad name
CC_FC_perc_noNeg = CC_FC_perc;
CC_FC_objL_noNeg = CC_FC_objL;
CC_FC_tree_noNeg = CC_FC_tree;

FC_sim_noNeg = FC_sim;
FC_sim_noNeg(find(FC_sim<0))=0;
for s=1:100
   CC_FC_noNeg(:,s) = clustering_coef_wu(FC_sim_noNeg(:,:,s)); 
end

d_FC_perc_noNeg = d_FC_perc;
d_FC_objL_noNeg = d_FC_objL;
d_FC_tree_noNeg = d_FC_tree;

d_FC_noNeg = d_FCsim_abs; %this one was mistakenly named

FC_sim_noNeg = FC_sim;
FC_sim_noNeg(find(FC_sim<0))=0;
for s=1:100
   d_FC_noNeg(s) = density_und(FC_sim_noNeg(:,:,s)); 
end

FC_sim_noNeg = FC_sim_abs; %this one was mistakenly named
FC_sim_perc_noNeg = FC_sim_perc;
FC_sim_objL_noNeg = FC_sim_objL;
FC_sim_tree_noNeg = FC_sim_tree;

NMI_sim_noNeg = NMI_sim;
NMI_sim_header_noNeg = NMI_sim_header;

%absolute value
load('Results_absValue.mat');

CC_FC_perc_abs = CC_FC_perc;
CC_FC_objL_abs = CC_FC_objL;
CC_FC_tree_abs = CC_FC_tree;

d_FC_perc_abs = d_FC_perc;
d_FC_objL_abs = d_FC_objL;
d_FC_tree_abs = d_FC_tree;
d_FC_abs = d_FCsim_abs;

FC_sim_perc_abs = FC_sim_perc;
FC_sim_objL_abs = FC_sim_objL;
FC_sim_tree_abs = FC_sim_tree;

NMI_sim_abs = NMI_sim;
NMI_sim_header_abs = NMI_sim_header;

clear CC_FC_objL CC_FC_perc CC_FC_tree d_FC_objL d_FC_perc d_FC_tree FC_sim_objL FC_sim_perc FC_sim_tree NMI_sim NMI_sim_header Modules_sim


%% Figure 3 Negative Weights

Edges_LFRw = Adj2lowerTriangleVector(LFR_orig);
numberEdges = size(Edges_LFRw,1);
Edges_FCsim = Tensor2EdgeMatrix(FC_sim);

e=1:numberEdges;%find(Edges_LFRw>0);
connected = Edges_FCsim(e,:);
min_connection_strength = min(connected,[],2);
figure;
scatter(Edges_LFRw(e),min_connection_strength,'filled','k','MarkerFaceAlpha',.5);
title('LFR edge weight versus population min FC_s_i_m weight');
xlabel('LFR weight before noise simulation');
ylabel('Population minimum correlation (R) of FC_s_i_m');
set(gca,'fontsize', 16) 


%% Figure 4  Objective Function Curve (for all or just lambda?)

%% Figure 4 NMI (bar)

NMI_DATA = [NMI_sim_abs(:,1),NMI_sim_abs(:,5), NMI_sim_abs(:,2),NMI_sim_noNeg(:,2),NMI_sim_abs(:,3),NMI_sim_abs(:,4),NMI_sim_abs(:,6),NMI_sim_noNeg(:,3),NMI_sim_noNeg(:,4),NMI_sim_noNeg(:,6)]
NMI_HEADER = {'FCsim','Statistical','Abs(FCsim)','NoNeg(FCsim)','Abs Percolation','Abs Objective(Lambda)','Abs Max Spanning Tree','No Neg Percolation','No Neg Objective(Lambda)','No Neg Max Spanning Tree'};

BarGraphs(NMI_DATA,...
    'Normalized Mutual Information: Binary FC_s_i_m versus LFR ground truth ',...
NMI_HEADER);

[h,p,ci,stats]=ttest(NMI_sim_abs(:,3),NMI_sim_abs(:,4))

[h,p,ci,stats]=ttest(NMI_sim_noNeg(:,3),NMI_sim_noNeg(:,4))

%% Figure 5 Densty (bar)
DENSITY_DATA = [d_FC_stat', d_FC_noNeg',d_FC_perc_abs',d_FC_objL_abs',d_FC_perc_noNeg',d_FC_objL_noNeg',d_FC_tree_noNeg']
D_HEADER = {'Statistical','NoNeg(FCsim)','Abs Percolation','Abs Objective(Lambda)','No Neg Percolation','No Neg Objective(Lambda)','No Neg Max Spanning Tree'};

BarGraphs(DENSITY_DATA - d_LFRb,...
    'Density Accuracy',...
D_HEADER);
ylabel('Density of FCsim - density of LFR')

[h,p,ci,stats]=ttest(d_FC_perc_abs,d_FC_objL_noNeg)


%% Figure 6 Clustering (bar)
CC_LFRb_rep = repmat(CC_LFRb,1,100,1);
CC_DATA = [mean(CC_FC_stat-CC_LFRb_rep,2), mean(CC_FC_noNeg-CC_LFRb_rep,2),mean(CC_FC_perc_abs-CC_LFRb_rep,2),mean(CC_FC_objL_abs-CC_LFRb_rep,2),mean(CC_FC_perc_noNeg-CC_LFRb_rep,2),mean(CC_FC_objL_noNeg-CC_LFRb_rep,2),mean(CC_FC_tree_noNeg-CC_LFRb_rep,2)];
CC_HEADER = {'Statistical','NoNeg(FCsim)','Abs Percolation','Abs Objective(Lambda)','No Neg Percolation','No Neg Objective(Lambda)','No Neg Max Spanning Tree'};

BarGraphs(CC_DATA,...
    'Clustering Coefficient Accuracy',...
CC_HEADER);
ylabel('CC of FCsim - CC of LFR')

[h,p,ci,stats]=ttest(mean(CC_FC_perc_abs-CC_LFRb_rep,2),mean(CC_FC_objL_abs-CC_LFRb_rep,2))

[h,p,ci,stats] = ttest(mean(CC_FC_stat-CC_LFRb_rep,2), mean(CC_FC_noNeg-CC_LFRb_rep,2))
