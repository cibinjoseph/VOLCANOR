dataMat = dlmread('Results/r01forceHist.txt');
dataMat(1,:) = [];
CT = dataMat(:,2);
indxStart = length(CT)-360;  % Last 5 revs
avgCT = mean(CT(indxStart:end));
disp(avgCT)
