function Nij = NMI_confusionEntry(ROIlistA,ROIlistB)
%derived from:
%Alexander-Bloch et al 2012 "The discovery of population
%difference in network community structure: New methods and applications to
%brain functional networks in schizophrenia" NeuroImage 59: 3889-3900

    %find the value of the entry in the confusion matrix based on the
    %number of ROIs each matrix has for the entry.
    Nij=0;

    LA=length(ROIlistA);
    LB=length(ROIlistB);
    if(LA<LB) %if A is the smaller, re-arrange
        ROIlist_temp=ROIlistA;
        ROIlistA=ROIlistB;
        ROIlistB=ROIlist_temp;
        clear ROIlist_temp
    end
    %re-define
    LA=length(ROIlistA); 
    LB=length(ROIlistB);

    diff=LA-LB;  %need to correct if not the same size

    padding=NaN(diff,1);
    ROIlistB=[ROIlistB;padding];

    for k=1:LA
       match=find(ROIlistB==ROIlistA(k));
       if(~isempty(match))
            Nij=Nij+1; 
       end  
    end

end

