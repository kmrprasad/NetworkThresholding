function [Sparsified] = threshold_pvalue(weightedNetworks,weightedNetworks_pvalue,alpha)
%Takes a weighted adjacency matrix (assumed to be a square matrix) and
%a edge-matched associated pValue network, and sets edge weight to zero in the weighted (correlation) network if the corresponding p-value in the pvalue network is below or equal to the
%alpha value 
%does not alter values above the threshold


%"weightedNetworks" is a N x N x M matrix, where N is the number of nodes and
%M is the number of networks

%weightedNetworks_pavlue is the same format, represneting the assocated
%p-values

%alpha is a scalar, the significant cutoff

%sparsNetwork is the same size as the original weighteNetworks, but with
%zeros where there were previously NaN or sub threshold values.
FCONp_mask = find(weightedNetworks_pvalue>alpha); %find the non-significant p-values

Sparsified=weightedNetworks;

Sparsified(FCONp_mask)=0;

end

