function [correctedMat] = avgDoubleCells(dubSessionInd,dataMat)

correctedMat = dataMat;

if min(size(dataMat)) > 1

for i = 1:length(dubSessionInd)
    
    correctedMat(dubSessionInd(i)-1,:) = mean(dataMat(dubSessionInd(i)-1:dubSessionInd(i),:));
    
end

correctedMat(dubSessionInd,:) = [];

else
    
    for i = 1:length(dubSessionInd)
        
        correctedMat(dubSessionInd(i)-1) = mean(dataMat(dubSessionInd(i)-1:dubSessionInd(i)));
        
    end
    
    correctedMat(dubSessionInd) = [];
    
end


