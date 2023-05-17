function Networks_Binary = threshold_percolation(Networks)
%PERCOLATION Summary of this function goes here
%Networks is the collection of adjacency matrices that you want to use, an
%R x R x N tensor where R = number Regions (nodes), and N = number Networks

%figure out if there is a third dimension
try
    numberNets = size(Networks,3);
catch
    numberNets = 1;
end    

%find dimensionality
numberNodes = size(Networks,1);
numberEdges=((numberNodes^2)-numberNodes)/2;

%preallocate the binary network output
Networks_Binary=zeros(numberNodes,numberNodes,numberNets);

%% Loop through all subjects
    for n = 1:numberNets

        % once again dealing with the possibility that third dimension is 1
        if (numberNets>1)
            Network_weighted=Networks(:,:,n);
        else
           Network_weighted=Networks(:,:);
        end

        % convert to vector, 
        EdgeVector_weighted = Tensor2EdgeMatrix(Network_weighted);
        
        %sort  ascending...
        [EdgeVector_weighted_ranked,rankOrder]=sort(EdgeVector_weighted,'Ascend');
        
        % create map of edgeIndex.
        [~,EdgeIndexMatrix] = EdgeVector2AdjacencyMatrix(EdgeVector_weighted_ranked,numberNodes);
             
        %initialize memory variable (loop will only find the percolation
        %threshold after its "too late"
        binary_previous=double(logical(Network_weighted)); %double(logical( means "binary network where all non-zero values are 1"

       
        %remove edges sequentially...
        for t = 1:numberEdges

          %set to previous run outcome ( set to initial values if t=1);
          binary_current=binary_previous;

          weakest_edge_index = rankOrder(t);
          
          %find the index of the new weakest edge, given by the t-th entry
          adjace_index_of_weakest_edge = find(EdgeIndexMatrix==weakest_edge_index);
          
          binary_current(adjace_index_of_weakest_edge)=0; %zero out new edges
          
          [~,temp_componentSizes] = get_components(binary_current); %get_components is a BCT function
          
          %figure out if all the nodes are in the component
          gcc_current=max(temp_componentSizes); %max is the largest component, should be == numberNodes if fully connected
        
          gcc_percent_current=gcc_current/numberNodes; %divide by size for percentage

          % now check if the GCC is still fully connected
          if(gcc_percent_current<1) %GCC NOT FULLY CONNECTED!!
                Networks_Binary(:,:,n)=binary_previous; %the previous iteration was the percolation threshold!!
                break; %stop executing edge loop
          else %if the GCC is still fully connected, i.e. %gcc=1
                 binary_previous=binary_current;
          end
        end %end of edge loop
    end %end of subject loop
    
    %clean up output
    Networks_Binary=squeeze(Networks_Binary); %removes third dimension if == 1
end %end of function
