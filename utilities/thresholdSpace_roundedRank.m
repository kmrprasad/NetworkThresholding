function [roundedRank_threshold_space] = thresholdSpace_roundedRank(weightedNetwork,precision)

    weightedEdgeVector = Adj2lowerTriangleVector(weightedNetwork); %get unique lower triangle off-diagonal elements
    weightedEdgeVector=round(weightedEdgeVector,precision); %round to the desired precision 3=thousandth
    roundedRank_threshold_space = unique(weightedEdgeVector); %remove duplicates
    roundedRank_threshold_space = sort(roundedRank_threshold_space,'Ascend'); %sort from lowest to highest

end

