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

gcc = 1;

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


Modules_sim=zeros(numberNodes,numberSim,6);
NMI_sim_header={'Weighted','Weighted No Negatives','Percolation','Objective Lambda',...
    'Statistical Threshold','Maximum Spanning Tree'};
NMI_sim = zeros(numberSim,6);

%%

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
    
      % Rounded Rank Threshold Space
    weightedEdgeVector = Adj2lowerTriangleVector(FC_sim_noNeg(:,:,h)); %get unique lower triangle off-diagonal elements
    weightedEdgeVector=round(weightedEdgeVector,3); %round to the nearest thousandth
    roundedRank_threshold_space = unique(weightedEdgeVector); %remove duplicates
    roundedRank_threshold_space(find(roundedRank_threshold_space<0))=[];
    roundedRank_threshold_space = sort(roundedRank_threshold_space,'Ascend'); %sort from lowest to highest

    
    [optima_threshold,optima_threshold_index,...
        optimum_F_values,header,M_all,F_all,M0_threshold_value,...
        M1_threshold_value,M0_threshold_indx,...
        M1_threshold_indx] = threshold_objective(temp,...
        roundedRank_threshold_space,gcc);
end  


%% plot
figure;
for g=1:length(header)
    subplot(3,3,g)
    plot(roundedRank_threshold_space,M_all(g,:),'k','linewidth',3);
    hold on;
     plot(roundedRank_threshold_space,F_all(g,:),'b','linewidth',3);

    xline(roundedRank_threshold_space(M0_threshold_indx),'r','linewidth',2);
    xline(roundedRank_threshold_space(M1_threshold_indx),'m','linewidth',2);

    xlim([roundedRank_threshold_space(M0_threshold_indx)-0.01,roundedRank_threshold_space(M1_threshold_indx)+0.01]);
    scatter(optima_threshold(g),F_all(g,optima_threshold_index(g)),[],'red','filled','d');
    legend({'M','F_M','M0','M1'});
    title(header(g));
    xlabel('threshold');
end
