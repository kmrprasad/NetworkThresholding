%LFR  Networks: https://github.com/carlonicolini/lfrwmx
%must download the source code from github
%must install CMAKE
%must have GCC compiler installed
%must compile source code to lfrwmx build
%probably will have make some symbolic links as indicated on the git page
%then must add build path to matlab, as below:
lfr_path='/data1/lfrwmx-master/build';
addpath(lfr_path);
%this will add the lfrw_mx function which is used below

%LFR parameters
lfr_nodes=363; %'N' number of nodes in network
lfr_degree_ave=24; %'k' average nodal degree
lfr_degree_max=200; %'maxk' maximum nodal degree
lfr_community_min_size = 12; %'minc' minimum community size
lfr_community_max_size = 32; %'maxc' maximum community size
lfr_community_mixing = 0.6; %'mut' community mixing coefficient [0,1]
lfr_weights_mixing = 0.2; %'muw' cweight mixing coefficient [0,1]
lfr_exp_degree = 2; %'t1' exponent of degree distribution
% lfr_exp_comm = .5; %'t2' exponenet for community size distribution
% lfr_beta = 1; %'beta' desired beta exponent
% lfr_CC = 5; %'C' desired clustering coefficient of network

%make lfr network
[network,community]=lfrw_mx('N',lfr_nodes,'k',lfr_degree_ave,'maxk',lfr_degree_max,...
     'mut',lfr_community_mixing,'muw',lfr_weights_mixing','t1',lfr_exp_degree);
    %'minc',lfr_community_min_size,'maxc',lfr_community_max_size,...
   

%plot results
figure;
imagesc(network);
colorbar;