function [binaryNetwork] = threshold_intensity(weightedNetwork,threshold)
%Takes a weighted adjacency matrix (assumed to be a square matrix) and
%a threshold value and sets edge weight to zero if below the threshold or one if above,
%resulting in a binray network
    weightedNetwork(isnan(weightedNetwork))=0; %clean up any NaN values if there are any
    weightedNetwork(weightedNetwork<=threshold)=0; %very important that this comes first
    weightedNetwork(weightedNetwork>threshold)=1; %see NOTE above: very important that this comes second
    binaryNetwork=weightedNetwork;
end
