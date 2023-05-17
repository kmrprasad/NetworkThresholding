%% Example_ObjectiveFunction
% for references and additional information, see README.md
% By: University of Pittsburgh CONCEPT LAB 2023

%% Configure the workspace and search paths
addpath('utilities');
Configure_LFR_sim; %adds folders and loads basic data

%% Get info on Loaded data:
%corMatOrig is the binary LFR simulation LFR network: n x n nodes
%corMatOrigw is the weighted LFR simulation LFR network: n x n nodes
%FC_sim are the noise-added replicates from the LFR simulation: 
   %    n x n nodes x r replicates
%m is the ground truth modularity of the LFR network: 1 x n nodes

numberNodes = size(corMatOrig,1);
numberReps = size(FC_sim,3);
 
%% Show objective function implementation on a single subject network:
GCC = 1; %assign the GCC value between (0,1] ...exclusive 0, inclusive 1
randomSubject = randi([1 100],1); %seelect a random subject
W = FC_sim(:,:,randomSubject); %specify the weighted network

disp(datetime);

disp('calculating threshold space, T');
%get the threshold space
T = thresholdSpace_roundedRank(W,3);

disp('finding binary network using objective function method');
% apply the objective function
[B,theta] = threshold_objective_LAMBDA(W,T,GCC); 

disp(datetime);

%% plot the result and compare the binary;
figure;

subplot(1,3,1);
imagesc(corMatOrig); colormap('bone');
colorbar; title('Original LFR Network');

subplot(1,3,2);
imagesc(W); colormap('jet');
colorbar; title('Noisy FC simulation on LFR Network');

subplot(1,3,3);
imagesc(B); colormap('bone');
colorbar; title('Binarized FC sim using Objective Function');




