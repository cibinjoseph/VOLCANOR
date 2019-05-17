clear; clf;
knp=dlmread('KnPdata.dat');
sim=dlmread('Results/r01forceHist.txt');
sim(1,:)=[];
plot(knp(:,1),knp(:,2),'ro');
hold on;
plot(sim(:,1)/16,sim(:,2),'b');
grid on;
