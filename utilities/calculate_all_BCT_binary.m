function [NetworkSummaryRow,header] = calculate_all_BCT_binary(network)

    if(nargin<1)

        header={'Density'; 'Average Nodal Degree';...
    'Characteristic Path Length'; 'Transitivity';
    'Average Clustering Coefficient'; 'Efficiency';...
    'Modularity';'Average Betweennes Centrality';'Assortativity'};
        NetworkSummaryRow=header;
        return;
    end

    r = assortativity_bin(network,0); 
    C = clustering_coef_bu(network);
    AveC=mean(C);
    try
    [Ci,Q] = modularity_und(network);
    catch
    Q=0;    
    end
    BC=betweenness_bin(network);
    AveBC=mean(BC);
    Dist=distance_bin(network);
    [lambda,efficiency,ecc,radius,diameter] = charpath(Dist,0,0); %no diagonal, no infinite distances
    d = degrees_und(network); 
    aved = mean(d); 
    [density, V, NumOfEdges] = density_und(network); 
    [T]=transitivity_bu(network);
    
    
    NetworkSummaryRow=[density,aved,lambda,T,AveC,efficiency,Q,AveBC,r];
    NetworkSummaryRow(isnan(NetworkSummaryRow))=0;
    NetworkSummaryRow=real(NetworkSummaryRow);
    
    
header={'Density'; 'Average Nodal Degree';...
    'Characteristic Path Length'; 'Transitivity';
    'Average Clustering Coefficient'; 'Efficiency';...
    'Modularity';'Average Betweennes Centrality';'Assortativity'};   
    

    
end
