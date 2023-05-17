function [BinaryNetworks,OptimalThresholds, M0_thresholds, M1_thresholds, M_values, F_values] = threshold_objective_LAMBDA(weightedNetworks,thresholds,gcc)

%uses an objective function to turn a threshold sweep of a graph measure, M, into a minimization
%problem over "F"

%assumes characteristic path length, "lambda" is the best choice of M

%calculates "F" between M0 and M1.
%M0 is defined at the lowest threshold where density is not 1
%M1 is defined at the highest threshold where the GCC >= the desired GCC

%INPUTS
%weightedNetworks" is the collection of adjacency matrices that you want to use, an
%R x R x N tensor where R = number Regions (nodes), and N = number weightedNetworks
%thresholds is a 1 x T length vector where T is the number of thresholds
%tested, for instance T=linspace(0,1,101);
%gcc is a scaler, expressed as a percent, the minimum giant connected
%component size used to set M1

try
    numberNets = size(weightedNetworks,3);
catch
    numberNets=1;
end
numberNodes= size(weightedNetworks,1);
numberThresholds = length(thresholds);

temp_BinaryNetworks=zeros(numberNodes,numberNodes,numberNets,numberThresholds);
BinaryNetworks=zeros(numberNodes,numberNodes,numberNets);

optimum_F_values = zeros(numberNets,1);
OptimalThresholds = zeros(numberNets,1);

optima_threshold_index = ones(numberNets,1);

Network_Density = zeros(numberNets,numberThresholds);
Network_GCC_perc = zeros(numberNets,numberThresholds);

% F and M
F_values=zeros(numberNets,numberThresholds);
M_values=zeros(numberNets,numberThresholds);

M0_values=zeros(numberNets,1);
M1_values=zeros(numberNets,1);
M0_thresholds=zeros(numberNets,1);
M1_thresholds=zeros(numberNets,1);

M0_threshold_indx=ones(numberNets,1);
M1_threshold_indx=ones(numberNets,1);

for n = 1:numberNets


  for t = 1:numberThresholds
      
    % Get the binary version of the subject leve weighted network at
    % threshold "i"
    Binary_=threshold_intensity(weightedNetworks(:,:,n),thresholds(t));
    temp_BinaryNetworks(:,:,n,t)=Binary_;
    
    Dist=distance_bin(Binary_);
    [lambda,~,~,~,~] = charpath(Dist,0,0); %no diagonal, no infinite distances
    Mi = lambda;
    
    M_values(n,t)=Mi; %add to master M

    % find the density of the network. used to define M0 threshold
    [density_, ~, ~] = density_und(Binary_); 
    Network_Density(n,t) =  density_;

    % find size of largest collection of nodes. used to define M1 threshold
    [~,temp_componentSizes] = get_components(Binary_); %get_components is a BCT function      
     maxCluster=max(temp_componentSizes);
    Network_GCC_perc(n,t)=maxCluster/numberNodes; %max(cc) will be the number of nodes in largest component.  divide by size for percentage
    
    %now check each of the above properties
    if(t>1) %if its not the first threshold (t=1 is always threshold=0)
        %define M_0 value
        if(Network_Density(n,t)<1 && Network_Density(n,t-1)==1) %if the threshold is the first to produce a graph with density less than 1            
            M0_threshold_indx(n)=t;
            M0_values(n)=M_values(n,t);
            M0_thresholds(n)=thresholds(t-1);
        end
        %define M_1 value
       if(Network_GCC_perc(n,t)<gcc && Network_GCC_perc(n,t-1)>=gcc)  %if this threshold is below the cutoff but the previous was not
            M1_threshold_indx(n)=t-1;   
            M1_values(n)=M_values(n,t-1);
            M1_thresholds(n)=thresholds(t);
       end    
    else %if its the first threshold tested
        %define M_0 values
        if(Network_Density(n,t)<1) %if the density is already less than 1 on the first threshold, then threshold zero is the M0 value              
            M0_threshold_indx(n)=t;  
            M0_values(n)=M_values(n,t);

        end
        %define M_1 value
        if(Network_GCC_perc(n,t)<gcc) %if the first value is not connected enough       
            M1_threshold_indx(n)=t;
            M1_values(n)=M_values(n,t);            
       end
    end
  end

%% Calculate Distance Measure, F, for all networks, M at thresholds i
        M0=M0_values(n);
        M1=M1_values(n);
        for t = 1:length(thresholds)
                %Get Mi from the master list M_values
                Mi=M_values(n,t);
                %Compare that Mi to the endpoints of the range to to get F, where: F(i)=(Mi-M0)^2 + (Mi-M1)^2.
                F_values(n,t)=(Mi-M0).^2 + (Mi-M1).^2; %this performs formula on entire vector of graph measures
        end
    

%% find optimum
              reasonable_range_values=F_values(n,M0_threshold_indx(n):M1_threshold_indx(n));
              reasonable_range_thresholds = thresholds(M0_threshold_indx(n):M1_threshold_indx(n));
              indx_range=M0_threshold_indx(n):M1_threshold_indx(n);
              
              if(length(reasonable_range_values)<=2) %if the range isn't a range, pick the lower
                 optima_threshold_index(n) = indx_range(M0_threshold_indx(n));
                 OptimalThresholds(n)=reasonable_range_thresholds(M0_threshold_indx(n));
                optimum_F_values(n)=reasonable_range_values(M0_threshold_indx(n));
                continue; %exit early
              end
                            
              [~,i_max] = max(reasonable_range_values);
              [~,i_min] = min(reasonable_range_values);
                       
              %special case 1: it's a minimization problem and not an endpoint
              if(i_min ~= 1 && i_min ~= length(reasonable_range_values))
                  optima_threshold_index(n) = indx_range(i_min);
                  OptimalThresholds(n)=reasonable_range_thresholds(i_min);
                  optimum_F_values(n)=reasonable_range_values(i_min);

              end
              
              %special case 2: it's a maximization problem and not an
              %endpoint
              if(i_max ~= 1 && i_max ~= length(reasonable_range_values))
                  optima_threshold_index(n) = indx_range(i_max);
                  OptimalThresholds(n)=reasonable_range_thresholds(i_max);
                  optimum_F_values(n)=reasonable_range_values(i_max);
              end
              
   BinaryNetworks(:,:,n)=temp_BinaryNetworks(:,:,n,optima_threshold_index(n));

end %end network loop

end %end function