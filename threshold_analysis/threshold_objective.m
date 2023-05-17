function [optima_threshold,optima_threshold_index,optimum_F_values,header,M_all,F_all,M0_threshold_value,M1_threshold_value,M0_threshold_indx,M1_threshold_indx] = threshold_objective(network_weighted,threshold_space,gcc)

if(nargin<3)
gcc_cutoff=.9; %was 0.7, then also 0.9
else
   gcc_cutoff = gcc; 
end
header={'Density'; 'Average Nodal Degree';...
    'Lambda'; 'Transitivity';
    'Average Clustering Coefficient'; 'Efficiency';...
    'Modularity';'Average Betweenness Centrality';'Assortativity'};


%Main outputs
optimum_F_values = zeros(9,1);
optima_threshold = zeros(9,1);
optima_threshold_index = zeros(9,1);

Network_Density = zeros(1,length(threshold_space));
Network_GCC_perc = zeros(1,length(threshold_space));

% F and M
F_all=zeros(9,length(threshold_space));
M_all=zeros(9,length(threshold_space));

M0_threshold_value=zeros(9,1);
M1_threshold_value=zeros(9,1);

M0_threshold_indx=1;
M1_threshold_indx=length(threshold_space);

  for t = 1:length(threshold_space)
      
    % Get the binary version of the subject leve weighted network at
    % threshold "i"
    Binary_=threshold_intensity(network_weighted(:,:),threshold_space(t));

    Mi=calculate_all_BCT_binary(Binary_)'; %calculate 9 summary graph measures from Binary Network
    M_all(:,t)=Mi; %add to master M

    % find the density of the network. used to define M0 threshold
    [density_, ~, ~] = density_und(Binary_); 
    Network_Density(t) =  density_;

    % find size of largest collection of nodes. used to define M1 threshold
    a = get_components(Binary_); %get_components is a BCT function
    c = arrayfun(@(x)length(find(a == x)), unique(a), 'Uniform', false); %gets the number of instances of each unique entry
    cc = cell2mat(c); %have to convert to matrix
    Network_GCC_perc(t)=max(cc)/size(Binary_,1); %max(cc) will be the number of nodes in largest component.  divide by size for percentage
    
    %now check each of the above properties
    if(t>1) %if its not the first threshold (t=1 is always threshold=0)
        %define M_0 value
        if(Network_Density(t)<1 && Network_Density(t-1)==1) %if the threshold is the first to produce a graph with density less than 1            
            M0_threshold_indx=t;
            M0_threshold_value(:)=M_all(:,t);

        end
        %define M_1 value
       if(Network_GCC_perc(t)<gcc_cutoff && Network_GCC_perc(t-1)>=gcc_cutoff)  %if this threshold is below the cutoff but the previous was not
            M1_threshold_indx=t-1;   
            M1_threshold_value(:)=M_all(:,t-1);


       end    
    else %if its the first threshold tested
        %define M_0 values
        if(Network_Density(t)<1) %if the density is already less than 1 on the first threshold, then threshold zero is the M0 value              
            M0_threshold_indx=t;  
            M0_threshold_value(:)=M_all(:,t);

        end
        %define M_1 value
        if(Network_GCC_perc(t)<gcc_cutoff) %if the first value is not connected enough       
            M1_threshold_indx=t;
            M1_threshold_value(:)=M_all(:,t);            
       end
    end
  end

%% Calculate Distance Measure, F, for all graph measure, M at thresholds i
    for g = 1:9   
        M0=M0_threshold_value(g);
        M1=M1_threshold_value(g);
        for t = 1:length(threshold_space)
                %Get Mi from the master list M_all
                Mi=M_all(g,t);
                %Compare that Mi to the endpoints of the range to to get F, where: F(i)=(Mi-M0)^2 + (Mi-M1)^2.
                F_all(g,t)=(Mi-M0).^2 + (Mi-M1).^2; %this performs formula on entire vector of graph measures
        end
    end

%% find optimum
    for g=1:9
              reasonable_range_values=F_all(g,M0_threshold_indx:M1_threshold_indx);
              reasonable_range_thresholds = threshold_space(M0_threshold_indx:M1_threshold_indx);
              indx_range=M0_threshold_indx:M1_threshold_indx;
              
              if(length(reasonable_range_values)<=2) %if the range isn't a range, pick the lower
                 optima_threshold_index(g) = indx_range(M0_threshold_indx);
                 optima_threshold(g)=reasonable_range_thresholds(M0_threshold_indx);
                optimum_F_values(g)=reasonable_range_values(M0_threshold_indx);
                continue; %exit early
              end
              

              
              [~,maxInd] = max(reasonable_range_thresholds);
              [~,minInd] = min(reasonable_range_thresholds);
              
              [~,i_max] = max(reasonable_range_values);
              [~,i_min] = min(reasonable_range_values);
              
              %set a default, then check special cases
              %default case: it's the mid-point
                            %get the middle inxed
%                midIndx = round((minInd+maxInd/2));
%                
%                optima_threshold_index(g) = indx_range(midIndx);
%                optima_threshold(g)=reasonable_range_thresholds(midIndx);
%                optimum_F_values(g)=reasonable_range_values(midIndx);
              
              %special case 1: it's a minimization problem and not an endpoint
              if(i_min ~= 1 && i_min ~= length(reasonable_range_values))
                  optima_threshold_index(g) = indx_range(i_min);
                  optima_threshold(g)=reasonable_range_thresholds(i_min);
                  optimum_F_values(g)=reasonable_range_values(i_min);

              end
              
              %special case 2: it's a maximization problem and not an
              %endpoint
              if(i_max ~= 1 && i_max ~= length(reasonable_range_values))
                  optima_threshold_index(g) = indx_range(i_max);
                  optima_threshold(g)=reasonable_range_thresholds(i_max);
                  optimum_F_values(g)=reasonable_range_values(i_max);
              end
              

    end


end