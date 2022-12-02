function NMI = NMI_calculation(ConfusionMatrix,numberROIs)
%derived from:
%Alexander-Bloch et al 2012 "The discovery of population
%difference in network community structure: New methods and applications to
%brain functional networks in schizophrenia" NeuroImage 59: 3889-3900
    C_A=size(ConfusionMatrix,1);
    C_B=size(ConfusionMatrix,2);

    N_RowSums=sum(ConfusionMatrix,1);
    N_ColSums=sum(ConfusionMatrix,2);
    N = numberROIs;
    
    C_numerator_sum=0;
    CA_denominator_sum=0;
    CB_denominator_sum=0;
    
    %numerator calculations
    for i=1:C_A
       for j=1:C_B 
            N_doti=N_ColSums(i);
            N_dotj=N_RowSums(j);
            
            C=ConfusionMatrix(i,j)*log((ConfusionMatrix(i,j)*N)/(N_doti*N_dotj));
            C(isnan(C))=0;
            
            C_numerator_sum=C_numerator_sum+C;
       end
    end
    %denominator calculations
    for i=1:C_A
         N_doti=N_ColSums(i);
         C=N_doti*log(N_doti/N);
         C(isnan(C))=0;
         
         CA_denominator_sum=CA_denominator_sum+C;
    end
    for j=1:C_B
         N_dotj=N_RowSums(j);
         C=N_dotj*log(N_dotj/N);
         C(isnan(C))=0;
         
         CB_denominator_sum=CB_denominator_sum+C;
    end   
    
    NMI = (-2*C_numerator_sum)/(CA_denominator_sum+CB_denominator_sum);
    
end

