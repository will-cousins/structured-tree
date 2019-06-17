%FUNCTION
%   getImpedanceLaplace - wrapper for impedanceLaplace

%INPUTS
%   s - value at which to evaluate laplace transform of impedance
%   rRoot - radius of vessel at root of tree
%   alpha,beta - vessel radius scaling parameters
%   lrr - length-to-radius ratio
%   rMin - minimum radius of tree
%   par - struct of parameters

%OUTPUST
%   Z - laplce transform of impedance evaluated at s

function [Z,table] = getImpedanceLaplaceMod(s,rRoot,alpha,beta,lrr,rMin,ZTermST,par)

maxGens = ceil(log(rMin / rRoot) / log(alpha(1))) + 1;
table = nan(maxGens);
[Z,table] = impedanceLaplaceMod(s,1,1,rRoot,alpha,beta,lrr,rMin,ZTermST,par,table);
