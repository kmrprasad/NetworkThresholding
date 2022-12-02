function [communityNMI, Ci1, Ci2, Q1, Q2] = NMI_NetMod_Main(network1,network2)
    %takes a pair of undirected (weighted or not, but will perform best on
    %binary networks), and oupputs 5 variables:
    %communityNMI: a value between 0 and 1 representing the overall agreement of the community (module) classifications of the two networks
    %   0 = no bits of information are shared (most different)
    %   1 =  all bits of information are shared (community classifications
    %   are identical)
    % Ci is the community assignment vector for each network (each node gets an integer)
    % Q is the modularity index value fore each network
   
    %this algorithm is implemented based on the following paper:
    
    %Alexander-Bloch et al 2012 "The discovery of population
    %difference in network community structure: New methods and applications to
    %brain functional networks in schizophrenia" NeuroImage 59: 3889-3900


    %the (binary or weighted) networks are inputs to the BCT modularity
    %funciton from the Brain Connectivity Toolbox.
    %this function gives two outputs: Ci, the community assignments (as
    %integers that assign groups to nodes), and Q, the modularity index
    %value, which is not needed but helpful to have.  Usually you want 
    %the networks to have similar Q values to each other
    try
    [Ci1,Q1]=modularity_und(network1); 
    [Ci2,Q2]=modularity_und(network2); 
    catch
       
       Ci1=ones(size(network1,1),1);
       Ci2=Ci1;
       Q1=0;
       Q2=0;
        
    end
    %The two community sturctures are then compared
    %This function relies on two other functions:
    %   NMI_confusionEntry.m ... this calculates the confusion matrix based
    %   on the community assignments
    %   NMI_calculation.m  ...this calculates the NMI based on the
    %   confusion matrix
    communityNMI = NMI_communityComparison(Ci1,Ci2);

end

