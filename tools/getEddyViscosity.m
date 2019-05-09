% Returns eddy viscosity coefficient value for provided rotor parameters
% Ref: Tauszig 1998

clc; clear;

% Input parameters 
radius = 1.143;  % in metres
omega_rpm = 130.89  % in rpm
density = 1.225  % in kg/m3
kinematicViscosity = 1.46*10^(-5)  % in m/s2
Nblades = 2
CT = 0.00460
a1 = 2*10^(-4);  % Ref. Bagai, Leishman, 2002


% Calculated parameters
omega = 2*pi*omega_rpm/60;
velTip = radius*omega;

% Compute net lift on a single blade
lift = CT*density*pi*radius^2*velTip^2/Nblades;

% Assuming triangular distribution along span for lift 
% with max lift per unit span(sectional lift) at 0.8R
Lmax=2*lift/radius;

% Finding freestream velocity at max L at 0.8R
velAtLmax = 0.8*radius*omega;

% Using Kutta-Joukowski to find gammaMax
gammaMax = Lmax/(density*velAtLmax)

Re_nu = gammaMax/kinematicViscosity

% Assuming strength of tip vortex to be gammaMax
eddyVisc = 1+a1*gammaMax/kinematicViscosity

