clc; clear; clf;
A = dlmread('forceDist.curve');
r = A(2:end,1);
fz = A(2:end,2);

radius = 1.143;
omega = 62.83;
rho = 1.2;
vtip = radius*omega;

rbar = r/radius;
ct = fz/(pi*rho*vtip*vtip*radius*radius)
plot(rbar,ct)
