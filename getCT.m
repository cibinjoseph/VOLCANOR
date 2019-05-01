clc; clear; clf;
A = dlmread('Results/r01forceHist.txt');
iter = A(2:end,1);
ct = A(2:end,2);

radius = 1.143;
omega = 130.90;
rho = 1.225;
vtip = radius*omega;
dt = 0.00138;

nrev = iter*dt*omega/(2*pi);
%rbar = r/radius;
%ct = fz/(pi*rho*vtip*vtip*radius*radius)
%plot(rbar,ct)
plot(nrev,ct)

