function getCL_PROWIM(avgFrom)
clc; clf;

if (length(avgFrom) == 1)
  avgFrom = [avgFrom(1) avgFrom(1)];
end

W=dlmread('Results/r01forceHist.txt');
W(1,:)=[];
subplot(2,1,1);
plot(1:size(W,1),W(:,2))
hold on;
line([avgFrom(1) avgFrom(1)],[min(W(:,2)) max(W(:,2))]);
hold off;


P=dlmread('Results/r02forceHist.txt');
P(1,:)=[];
subplot(2,1,2);
plot(1:size(P,1),P(:,3))
hold on;
line([avgFrom(2) avgFrom(2)],[min(P(:,3)) max(P(:,3))]);

format long;
avgW=mean(W(avgFrom(1):end,2))
avgP=mean(P(avgFrom(2):end,3))
return;

