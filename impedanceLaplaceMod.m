%FUNCTION
%   impedanceLaplace - computes Z(s), laplace transform of impedance, using 
%       ST algorithm (see Cousins, Gremaud, 2013) structured tree paper for 
%       more details. no tiers here. difference between this and
%       impedanceGremLaplce is that this code includes radius dependent
%       viscosity and is implemented with a terminal pressure

%INPUTS
%   s - value at which to evaluate laplace transform of impedance
%   nAlpha, nBeta - identifies location of vessel within tree. should be set to 1
%       when calling function
%   alphaST, betaST - scaling parameters
%   rRoot - radius of root vessel of tree
%   rMinST - minimum radius
%   zTerm - terminal impedance
%   k1,k2,k3 - parameters for pressure/area constitutive law
%   lrrST - length-to-radius ratio
%   par - struct of parameters (needs rho, mu, k(1),k(2),k(3), gamma)
%   table - table of impedance values in all vessels of tree. should be set
%       to a matrix of NaN when calling function

%OUTPUTS
%   Z0 - impedance at root of tree
%   table - table of impedances for all vessels in tree

function [Z0,table] = impedanceLaplaceMod(s,nAlpha,nBeta,rRoot,alphaST,betaST,lrrST,rMinST,ZTermST,par,table)

r = rRoot * alphaST ^ (nAlpha - 1) * betaST ^ (nBeta - 1);
diam = 2 * r;
diam = 1e4 * diam; %convert cm to microns
etaStar = 6*exp(-0.085*diam)+3.2-2.44*exp(-0.06*diam.^.645);
mFactor = par.mu / 3.2;
muStar = mFactor*(1 + (etaStar-1).*(diam./(diam-1.1)).^2).*(diam./(diam-1.1)).^2;

A0 = pi * r ^ 2;
L = lrrST * r;
Ehoverr = par.k(1) * exp(par.k(2) * r) + par.k(3);
C = 3 * A0 / (2 * Ehoverr);
delta = 2 * (par.gamma + 2) * muStar / (par.rho * r ^ 2);
d = sqrt(A0 / (s * (s + delta) *  par.rho * C));

if r < rMinST
    ZL = ZTermST;
else
    if isnan(table(nAlpha+1,nBeta))
        [ZD1,table] = impedanceLaplaceMod(s,nAlpha+1,nBeta,rRoot,alphaST,betaST,lrrST,rMinST,ZTermST,par,table);
    else
        ZD1 = table(nAlpha+1,nBeta);
    end
    
    if isnan(table(nAlpha,nBeta+1))
        [ZD2,table] = impedanceLaplaceMod(s,nAlpha,nBeta+1,rRoot,alphaST,betaST,lrrST,rMinST,ZTermST,par,table);
    else
        ZD2 = table(nAlpha,nBeta+1);
    end
    ZL = (ZD1 * ZD2) / (ZD1 + ZD2);
end

if s == 0
    Z0 = ZL + 2 * (par.gamma + 2) * muStar * lrrST / (pi * r ^ 3);
else    
    num = ZL + tanh(L / d) / (s * d * C);
    denom = s * d * C * ZL * tanh(L / d) + 1;
    Z0 = num / denom;
end

table(nAlpha,nBeta) = Z0;