% creates figure 5.1 in A NEW PHYSIOLOGICAL BOUNDARY CONDITION FOR
% HEMODYNAMICS
r = 0.11;
Nt = 40;
dt = 0.025;
alpha = 0.91;
beta = 0.58;
lrr = 50;
rMin = 0.0083;
ZTerm = 0;

par = getParams;

z = lubichCoefMod(Nt,dt,r,alpha,beta,lrr,rMin,ZTerm,par);

figure; clf
plot(z,'.-')
axis tight,set(gca,'YLim',[-1000,12500])