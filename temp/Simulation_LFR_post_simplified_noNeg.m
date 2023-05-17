%% To be run after the R simulation
%%
% allocate results

Configure_LFR_sim;


%% load in LFR  data
Analysisfilepath='/data1/LFR_MRI_sim/scripts';

%real LFR network
LFR_data = sprintf('%s/PowerLawNetMat.mat',Analysisfilepath);
load(LFR_data);
modules_groundtruth = m'; %transpose
node_number = length(modules_groundtruth);

numberNodes=size(corMatOrig,1);
numberEdges = (numberNodes^2)-numberNodes;


%% noisy  FCN
time_start = datetime;

GCC = 1;

% read in R outputs
simulation_data_folder = sprintf('%s/Data',Analysisfilepath);
simulation_subject_data_dir = sprintf('%s/Simulation_100Subj_SNR35_Try1_Mut2',simulation_data_folder);

simList=dir(simulation_subject_data_dir);
simList(1:2)=[];




% allocation
numberSim=size(simList,1);

LFR_orig = corMatOrigw;
LFR_sim_numberNeg=zeros(numberSim,2);

threshold_value_Perc=zeros(numberSim,1);
threshold_value_ObjL=zeros(numberSim,1);

FC_sim = zeros(numberNodes,numberNodes,numberSim);
FC_sim_noNeg = zeros(numberNodes,numberNodes,numberSim);
FC_sim_perc = zeros(numberNodes,numberNodes,numberSim);
FC_sim_objL = zeros(numberNodes,numberNodes,numberSim);
FC_sim_stat = zeros(numberNodes,numberNodes,numberSim);
FC_sim_tree = zeros(numberNodes,numberNodes,numberSim);


CC_LFRb = clustering_coef_bu(corMatOrig);
CC_FCsim = zeros(numberNodes,numberSim);
CC_FCsim_noNeg = zeros(numberNodes,numberSim);
CC_FC_perc = zeros(numberNodes,numberSim);
CC_FC_objL = zeros(numberNodes,numberSim);
CC_FC_stat = zeros(numberNodes,numberSim);
CC_FC_tree = zeros(numberNodes,numberSim);


d_LFRb = density_und(corMatOrig);
d_FCsim = zeros(1,numberSim);
d_FCsim_noNeg = zeros(1,numberSim);
d_FC_perc = zeros(1,numberSim);
d_FC_objL = zeros(1,numberSim);
d_FC_stat = zeros(1,numberSim);
d_FC_tree = zeros(1,numberSim);


Modules_sim=zeros(numberNodes,numberSim,6);
NMI_sim_header={'Weighted','Weighted No Negatives','Percolation','Objective Lambda',...
    'Statistical Threshold','Maximum Spanning Tree'};
NMI_sim = zeros(numberSim,6);

for h = 1:numberSim
    %load data
    simName=simList(h).name; %get the subject ID
    matfile = sprintf('%s/%s',simulation_subject_data_dir,simName);
    load(matfile);
    
    disp('loaded replicate');
    disp(h);
    disp(datetime);
    
    temp_corr = RealCorMat - diag(diag(RealCorMat)); %zero out the diagonal
   
     FC_sim(:,:,h)=temp_corr;
    
    temp=temp_corr;
    temp(find(temp<0))=0;
       
    FC_sim_noNeg(:,:,h)=temp;
    
    %weighted
    disp('weighted');
    modules_estimate_1 = modularity_und(FC_sim(:,:,h));
    Modules_sim(:,h,1)= modules_estimate_1;
    NMI_sim(h,1) = NMI_communityComparison(modules_groundtruth,modules_estimate_1);
    LFR_sim_numberNeg(h,1)=numel(find(FC_sim(:,:,h)<0));
    LFR_sim_numberNeg(h,2)=LFR_sim_numberNeg(h,1)/numberEdges;
    disp(LFR_sim_numberNeg(h,2));
    
    %absolute value
    disp('abs of weightede');
    modules_estimate_2 = modularity_und(FC_sim_noNeg(:,:,h));
    Modules_sim(:,h,2)= modules_estimate_2;
    NMI_sim(h,2) = NMI_communityComparison(modules_groundtruth,modules_estimate_2);
    
    
    % Rounded Rank Threshold Space
    weightedEdgeVector = Adj2lowerTriangleVector(FC_sim_noNeg(:,:,h)); %get unique lower triangle off-diagonal elements
    weightedEdgeVector=round(weightedEdgeVector,3); %round to the nearest thousandth
    roundedRank_threshold_space = unique(weightedEdgeVector); %remove duplicates
    roundedRank_threshold_space(find(roundedRank_threshold_space<0))=[];
    roundedRank_threshold_space = sort(roundedRank_threshold_space,'Ascend'); %sort from lowest to highest
%   fixed_threshold_space = linspace(-1,1,81);
%   fixed_threshold_space = linspace(0,1,41);

    %percolation threshold
    disp('percolation threshold')
    %FC_sim_perc(:,:,h) = threshold_percolation(FC_sim(:,:,h)); 
    [FC_sim_perc(:,:,h),threshold_value_Perc(h)] = threshold_percolation_depricated(FC_sim_noNeg(:,:,h),roundedRank_threshold_space);
    module_estimate_3 =  modularity_und(FC_sim_perc(:,:,h));
    Modules_sim(:,h,3) = module_estimate_3;
    NMI_sim(h,3) = NMI_communityComparison(modules_groundtruth,module_estimate_3);
 

    % objective threshold lambda
    disp('objective threshold');
    [FC_sim_objL(:,:,h),threshold_value_ObjL(h)] = threshold_objective_LAMBDA(FC_sim_noNeg(:,:,h),roundedRank_threshold_space,GCC);
    module_estimate_4 = modularity_und(FC_sim_objL(:,:,h));
    Modules_sim(:,h,4)= module_estimate_4;
    NMI_sim(h,4) = NMI_communityComparison(modules_groundtruth,module_estimate_4);

    disp('statistical threshold');
    %statistical threshold
    p_value_matrix=zeros(node_number,node_number);
    fn = fieldnames(RealTS);
    for k=1:numel(fn)
        if( isnumeric(RealTS.(fn{k})) )
            for c=1:numel(fn)
                if( isnumeric(RealTS.(fn{c})) )
                [~,p]=corrcoef(RealTS.(fn{k}),RealTS.(fn{c}));
                                p_value_matrix(k,c) = p(1,2);
                end
            end            
        end
    end
    FC_sim_stat(:,:,h)=threshold_pvalue(RealCorMat,p_value_matrix,0.05);
    modules_estimate_5 = modularity_und( FC_sim_stat(:,:,h));
    Modules_sim(:,h,5)= modules_estimate_5;
    NMI_sim(h,5) = NMI_communityComparison(modules_groundtruth,modules_estimate_5);

    
    
        % tree threshold lambda
    disp('max tree');
    FC_sim_tree(:,:,h) = UndirectedMaximumSpanningTree(FC_sim_noNeg(:,:,h));
    module_estimate_6 = modularity_und(FC_sim_tree(:,:,h));
    Modules_sim(:,h,6)= module_estimate_6;
    NMI_sim(h,6) = NMI_communityComparison(modules_groundtruth,module_estimate_6);
    
    %basset sigma range
%     [sigma, sigma_range] = BassetSmallWorldness(RealCorMat,threshold_space);
%     
%     sigma_low = threshold_space(sigma_range(1));
%     sigma_high = threshold_space(sigma_range(2));
%  
%     modules_estimate_4 = modularity_und(threshold_intensity(RealCorMat,sigma_low));
%     modules_estimate_5 = modularity_und(threshold_intensity(RealCorMat,sigma_high));
       
    %objective function
%     [optima_threshold,optima_threshold_index,optimum_F_values,header,M_all,F_all,M0,M1] = threshold_objective(RealCorMat,threshold_space);
%     %add to plot
%     plot(threshold_space,F_all(8,:),'r');
%     scatter(optima_threshold(8),optimum_F_values(8),'k');
%     plot(threshold_space,M_all(8,:),'b');
%     scatter(optima_threshold(8),M_all(1,optima_threshold_index(8)),'k');
%     
%     modules_estimate_6 = modularity_und(threshold_intensity(RealCorMat,optima_threshold(1))); 
%     modules_estimate_7 = modularity_und(threshold_intensity(RealCorMat,optima_threshold(2))); 
%     modules_estimate_8 = modularity_und(threshold_intensity(RealCorMat,optima_threshold(3))); 
%     modules_estimate_9 = modularity_und(threshold_intensity(RealCorMat,optima_threshold(4))); 
%     modules_estimate_10 = modularity_und(threshold_intensity(RealCorMat,optima_threshold(5))); 
%     modules_estimate_11 = modularity_und(threshold_intensity(RealCorMat,optima_threshold(6))); 
%     modules_estimate_12 = modularity_und(threshold_intensity(RealCorMat,optima_threshold(7))); 
%     modules_estimate_13 = modularity_und(threshold_intensity(RealCorMat,optima_threshold(8))); 
%     modules_estimate_14 = modularity_und(threshold_intensity(RealCorMat,optima_threshold(9))); 
%     


CC_FCsim(:,h) = clustering_coef_wu(FC_sim(:,:,h));
CC_FCsim_noNeg(:,h) =  clustering_coef_wu( FC_sim(:,:,h));
CC_FC_perc(:,h) = clustering_coef_bu( FC_sim_perc(:,:,h));
CC_FC_objL(:,h) = clustering_coef_bu( FC_sim_objL(:,:,h));
CC_FC_stat(:,h) = clustering_coef_bu( FC_sim_stat(:,:,h));
CC_FC_tree(:,h) = clustering_coef_bu( FC_sim_tree(:,:,h));



d_FCsim(:,h) = density_und(FC_sim(:,:,h));
d_FCsim_noNeg(:,h) =  density_und( FC_sim_noNeg(:,:,h));
d_FC_perc(:,h) = density_und( FC_sim_perc(:,:,h));
d_FC_objL(:,h) = density_und( FC_sim_objL(:,:,h));
d_FC_stat(:,h) = density_und( FC_sim_stat(:,:,h));
d_FC_tree(:,h) = density_und( FC_sim_tree(:,:,h));

end

time_complete = datetime;


BarGraphs(NMI_sim,'NMI',NMI_sim_header);

density=[d_FCsim',d_FCsim_noNeg',d_FC_perc',d_FC_objL',d_FC_stat',d_FC_tree'];
density_difference=d_LFRb-density;


%% Stats

%NMI
[tstatsNMI_weighted_v_perc,tstats_header]=stats_TtestAllStats(NMI_sim(:,2),NMI_sim(:,3))

[tstatsNMI_weighted_v_obj,~]=stats_TtestAllStats(NMI_sim(:,2),NMI_sim(:,4))

[tstatsNMI_perc_v_obj,~]=stats_TtestAllStats(NMI_sim(:,3),NMI_sim(:,4)) %percolation versus objective

%clustering
CC_distance_perc = zeros(numberNodes,numberSim);
CC_distance_objL = zeros(numberNodes,numberSim);

for h=1:numberSim
    CC_distance_perc(:,h)=CC_LFRb-CC_FC_perc(:,h);
    CC_distance_objL(:,h)=CC_LFRb-CC_FC_objL(:,h);
end

tstatsCC_distance_perc_v_objL=zeros(numberNodes,11);
for n=1:numberNodes
    tstatsCC_distance_perc_v_objL(n,:)=stats_TtestAllStats(CC_FC_perc(n,:),CC_FC_objL(n,:));
end



%% density...
figure
plot(density_difference,'LineWidth',5);
lgd=legend(NMI_sim_header);
xlabel('replicate'); ylabel('LFR density - binary density');
title('Accuracy of Density Value after Thresholding');
lgd.FontSize = 16;
ax = gca; 
ax.FontSize = 16; 

%% clustering...

figure; hold on;
plot(CC_LFRb - mean(CC_FCsim,2),'Color',[.4 .4 .4],'LineWidth',3)
plot(CC_LFRb - mean(CC_FCsim_abs,2),'Color',[.8 .8 .8],'LineWidth',2) % r for red
plot(CC_LFRb - mean(CC_FC_perc,2),'Color','r','LineWidth',3) % r for red
plot(CC_LFRb - mean(CC_FC_objL,2),'Color','b','LineWidth',2) % r for red
plot(CC_LFRb - mean(CC_FC_stat,2),'Color',[.6 .2 .6],'LineWidth',2) % r for red
plot(CC_LFRb - mean(CC_FC_tree,2),'Color',[.2 .8 .2],'LineWidth',3) % r for red
lgd=legend(NMI_sim_header);
xlabel('node'); ylabel('group average LFR CC - binary CC');
title('Accuracy of CC Value after Thresholding');
lgd.FontSize = 16;
ax = gca; 
ax.FontSize = 16; 

%% negative edges and zero weights...

Edges_LFRw = Adj2lowerTriangleVector(LFR_orig);
numberEdges = size(Edges_LFRw,1);
Edges_FCsim = Tensor2EdgeMatrix(FC_sim);

e=1:numberEdges%find(Edges_LFRw>0);
connected = Edges_FCsim(e,:);
min_connection_strength = min(connected,[],2);
figure;
scatter(Edges_LFRw(e),min_connection_strength,'filled','k');
title('Non-zero LFR edges versus population min FCsim');
xlabel('LFR weight before noise');
ylabel('Population Minimum Correlation (R) of FCsim (noisy) edge');

e_z=find(Edges_LFRw==0);
connected_z = Edges_FCsim(e_z,:);
min_connection_strength_zero = min(connected_z,[],2);
figure;
histogram(min_connection_strength);
hold on;
histogram(min_connection_strength_zero);
legend('Non-zero LFRw', 'Unconnected LFRw');
title('Edge weight histograms of population minimum FCsim by LFR connectivity')
xlabel('population minimum FC-sim weight (R)');
ylabel('number of edges with weight');

numel(find(mean(connected,2)<0))