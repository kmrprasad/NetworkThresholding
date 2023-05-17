function [boolout,percent_gcc,size_gcc,components] = isFullyConnected(network_adjacencyMatrix)

        numberNodes=size(network_adjacencyMatrix,2);

        [components,temp_componentSizes] = get_components(network_adjacencyMatrix); %get_components is a BCT function
          
        %figure out if all the nodes are in the component
        size_gcc=max(temp_componentSizes); %max is the largest component, should be == numberNodes if fully connected
        
        percent_gcc=size_gcc/numberNodes; %divide by size for percentage
        
        if(percent_gcc==1)
           boolout=true;
        else
            boolout=false;
        end

end

