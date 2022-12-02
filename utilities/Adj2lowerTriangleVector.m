function [vectorizedLowerTriangle] = Adj2lowerTriangleVector(AdjacencyMatrix)
%ADJ2UPPERTRIANGLEVECTOR Summary of this function goes here
%   the adjacency matrix is a square matrix (not tensor)
    numberROIs = size(AdjacencyMatrix,1);
    %use this to remove duplicates and diagonals
    A = ones(numberROIs);
    B = tril(A,-1); %extract only the elements below the main diagonal
    C=reshape(B,numel(B),1);
    zeroElements=find(~C); %find the zeros

    subjectNetwork=tril(AdjacencyMatrix);
    vectorizedLowerTriangle=reshape(subjectNetwork,numel(subjectNetwork),1);
    vectorizedLowerTriangle(zeroElements)=[];%remove the diagonal and lower triangle
     
   %index test
%      a = [1:81]';
%      ra = repmat(a,1,81);
%      for i=1:81
%        ra(:,i) = ra(:,i) + (81*(i-1)) ;
%      
%      end
end

