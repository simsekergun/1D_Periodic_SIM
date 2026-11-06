% RT vs Theta
clear;
% if you don't want parallel computing, comment off the lines 5-12 and 
% use "for" instead of "parfor" 
parpool_no = 8;
p = gcp('nocreate');
if isempty(p)
    parpool(parpool_no);
elseif p.NumWorkers ~= parpool_no
    delete(gcp('nocreate'));
    parpool(parpool_no);
end

a = 120;        % radius (nm)
geom = struct('type','circle','a',a);    
d    = 400;                                         % period (nm)
lams = 550:750;                 % wavelengths
theta_i = 55*pi/180;            % incidence angles
Nsim = 128;

R = zeros(length(lams),1);
T = R;
A = R;
res = R;
mR = R;
mT = R;
mA = R;

kt = 1;
parfor k = 1:length(lams)
    ep = eps_gold(lams(k));
    m = mpm1d(a, lams(k), ep, d, struct('pol','TE','thetai',theta_i));
    o = sim1dperiodic(geom, lams(k), ep, Nsim, d, struct('pol','TE','thetai',theta_i));
    R(k) = o.R;
    T(k) = o.T;
    A(k) = o.A;
    res(k) = o.energy_residual;
    mR(k) = m.R;
    mT(k) = m.T;
    mA(k) = m.A;
end

fig = figure(7); clf;
plot(lams, T, lams, mT,'--',lams, R, lams, mR,'--');
ylabel('R and T');
legend('T (SIM)', 'T (MPM)','R (SIM)', 'R (MPM)','Location','Best')
xlabel('Wavelength (nm)');
grid on
