clc; clear; 
A=dlmread('Results/r01forceHist.txt');
A(1,:)=[];
CL=mean(A(300:700,2))
