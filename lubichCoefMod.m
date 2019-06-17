%FUNCTION
%   lubichCoefMod - computes convolution quadrature coefficients for the
%       impedance. the "Mod" means this algorithm is very slightly modified
%       from that presented in Cousins et. al. 2013.  Only difference here
%       is that we automatically zero weights corresponding to t > 1,
%       rather than numerically compute them (weights decay to zero very
%       quickly so we lose little accuracy by doing this)

%INPUTS
%   Nt - max number of time steps in simulation
%   dt - time step size
%   rRoot - root radius of structured tree
%   alpha, beta - scaling coefficients for tree
%   lrr - length to radius ratio
%   rMin - minimum radius of structured tree
%   par - blood parameters (viscosity, density, etc)

%OUTPUTS
%   z - coefficients for impedance convolution quadrature

function z = lubichCoefMod(Nt,dt,rRoot,alpha,beta,lrr,rMin,ZTerm,par)

NtMod = min(Nt,round(1/dt));
NInt = 2 * NtMod;

eps = 1e-12;
rho = eps ^ (1 / (NInt));

theta = linspace(0,2*pi,NInt+1).';
xi = rho * exp(1i * theta);
delta = 3/2 - 2 * xi + xi.^2 /2;

z = zeros(Nt+1,1);
ZEval = zeros(size(delta));
for ii = 1:length(delta)
    ZEval(ii) = getImpedanceLaplaceMod(delta(ii) / dt,rRoot,alpha,beta,lrr,rMin,ZTerm,par);
end

trapWeights = [1;2*ones(NInt-1,1);1] / (2 * NInt);
for ii = 0:NtMod
    z(ii+1) = sum(trapWeights .* ZEval .* exp(-1i * ii * theta));
    z(ii+1) = z(ii+1) / (rho ^ ii);
end
z = real(z);