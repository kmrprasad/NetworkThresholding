%% configure

cd ../
Configure_LFR_sim;
cd threshold_analysis

%% preallocate

numberReplicates = size(FC_sim,3);

BinaryNetworks = cell(numberReplicates,1);
OptimalThresholds = zeros(9,numberReplicates);
M0_thresholds = zeros(9,numberReplicates);
M1_thresholds= zeros(9,numberReplicates);
M_values = cell(numberReplicates,1);
F_values = cell(numberReplicates,1);

NMI_all = zeros(numberReplicates,9);

%% evaluate

gcc=1;

for r=1:numberReplicates

    disp(r);
    
    weightedNetworks = FC_sim(:,:,r);
    weightedNetworks = abs(weightedNetworks);
    threshold_space = thresholdSpace_roundedRank(weightedNetworks,3);

    %[BinaryNetworks,OptimalThresholds, M0_thresholds, M1_thresholds, M_values, F_values] = threshold_objective_LAMBDA(weightedNetworks,threshold_space,gcc);

    [header,BinaryNetworks{r},OptimalThresholds(:,r),M0_thresholds(:,r),M1_thresholds(:,r),M_values{r},F_values{r}] = threshold_objective(weightedNetworks,threshold_space,gcc);
       
    BinaryNetwork_temp = BinaryNetworks{r};
    for g=1:9
        modules_calculated = modularity_und(BinaryNetwork_temp(:,:,g));
        NMI_all(r,g) = NMI_communityComparison(modules_groundTruth,modules_calculated); 
    end
end

%% collate
ix_ofInterest = [1,2,3,4,5,6];

%% plot
BarGraphs(NMI_all(:,ix_ofInterest),'',header(ix_ofInterest));

[h,p,ci,stats]=ttest(NMI_all(:,3),NMI_all(:,4));
