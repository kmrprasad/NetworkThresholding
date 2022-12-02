function [NMI,ConfussionMatrix] = NMI_communityComparison(Ci1,Ci2)
%derived from:
%Alexander-Bloch et al 2012 "The discovery of population
%difference in network community structure: New methods and applications to
%brain functional networks in schizophrenia" NeuroImage 59: 3889-3900
    numberROIs=size(Ci1,1);
    numberMods_1=max(Ci1);
    numberMods_2=max(Ci2);
    ConfussionMatrix=zeros(numberMods_1,numberMods_2);

    %now loop through the moduels to create the confusion matrix
    for i=1:numberMods_1
        for j=1:numberMods_2

            Ci=find(Ci1==i);
            Cj=find(Ci2==j);

            Nij=NMI_confusionEntry(Ci,Cj);
            ConfussionMatrix(i,j)=Nij;

        end
    end %end confusion loop

     NMI = NMI_calculation(ConfussionMatrix,numberROIs);
            
end


