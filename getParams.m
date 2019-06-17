%FUNCTION
%   getParams - outputs physical parameters used in simulation

%INPUTS
%   none

%OUTPUTS
%   par - struct of parameters

function par = getParams

par.mmHgtoBa = 1.33322e3;
par.mu = 0.048;
par.gamma = 2;
par.gRat = (par.gamma + 2) / (par.gamma + 1);
k1 = 2e7; 
k2 = -22.53;
k3 = 8.65e5;
par.k = [k1;k2;k3];
par.rho = 1.055;
par.p0 = 100 * par.mmHgtoBa;
