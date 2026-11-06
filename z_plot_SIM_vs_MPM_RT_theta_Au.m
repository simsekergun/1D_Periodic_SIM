% RT vs Theta
clear;
set(0, 'DefaultFigureCreateFcn', @(h,~) set(h, 'Renderer', 'painters'))
% if you don't want parallel computing, comment off the lines 6-13 and 
% use "for" instead of "parfor" 
parpool_no = 8;
p = gcp('nocreate');
if isempty(p)
    parpool(parpool_no);
elseif p.NumWorkers ~= parpool_no
    delete(gcp('nocreate'));
    parpool(parpool_no);
end

a = 130;            % radius
geom = struct('type','circle','a',a);    % nm
d    = 400;                                         % period (nm)
lams = 650;         % wavelenth (nm)
thetas = 0:80;      % incidence angles
Nsim = 128;     % SIM mesh sampling number

R = zeros(length(thetas),1);
T = R;
A = R;
res = R;
mR = R;
mT = R;
mA = R;

Rtm = zeros(length(thetas),1);
Ttm = R;
Atm = R;
restm = R;
mRtm = R;
mTtm = R;
mAtm = R;

kt = 1;
parfor k = 1:length(thetas)
    ep = eps_gold(lams);
    m = mpm1d(a, lams, ep, d, struct('pol','TE','thetai',thetas(k)*pi/180));
    o = sim1dperiodic(geom, lams, ep, Nsim, d, struct('pol','TE','thetai',thetas(k)*pi/180));
    mTM = mpm1d(a, lams, ep, d, struct('pol','TM','thetai',thetas(k)*pi/180));
    oTM = sim1dperiodic(geom, lams, ep, Nsim, d, struct('pol','TM','thetai',thetas(k)*pi/180));
    R(k) = o.R;
    T(k) = o.T;
    A(k) = o.A;
    res(k) = o.energy_residual;
    mR(k) = m.R;
    mT(k) = m.T;
    mA(k) = m.A;
    %
    Rtm(k) = oTM.R;
    Ttm(k) = oTM.T;
    Atm(k) = oTM.A;
    restm(k) = oTM.energy_residual;
    mRtm(k) = mTM.R;
    mTtm(k) = mTM.T;
    mAtm(k) = mTM.A;
end

fig = figure(7); clf;
plot(thetas, T, thetas, mT,'s',thetas, R, thetas, mR,'o');
ylabel('R and T for TE_z');
legend('T (SIM)', 'T (MPM)','R (SIM)', 'R (MPM)','Location','NorthEast')
xlabel('Incidence Angle (degrees)');
grid on

fig = figure(8); clf;
plot(thetas, Ttm, thetas, mTtm,'s',thetas, Rtm, thetas, mRtm,'o');
ylabel('R and T for TM_z');
legend('T (SIM)', 'T (MPM)','R (SIM)', 'R (MPM)','Location','East')
xlabel('Incidence Angle (degrees)');
grid on

fig = figure(9); clf;
fig.Position = [100 100 500 700];
subplot(211);
plot(thetas, T, thetas(1:3:end), mT(1:3:end),'x',thetas, R, thetas(1:3:end), mR(1:3:end),'o');
ylabel('{\it{R}} and {\it{T}}');
text(5, 0.5, ' TE_z ','EdgeColor', 'black')
legend('{\it{T}} (SIM)', '{\it{T}} (MPM)','{\it{R}} (SIM)', '{\it{R}} (MPM)','Location','NorthEast')
xlabel('Incidence Angle (degrees)');
text(-12, 1.05, '(a)', 'FontSize',24);
grid on
%
subplot(212);
plot(thetas, Ttm, thetas(1:3:end), mTtm(1:3:end),'x',thetas, Rtm, thetas(1:3:end), mRtm(1:3:end),'o');
ylabel('{\it{R}} and {\it{T}}');
text(5, 0.5, ' TM_z ','EdgeColor', 'black')
text(-12, 1.05, '(b)', 'FontSize',24);
legend('{\it{T}} (SIM)', '{\it{T}} (MPM)','{\it{R}} (SIM)', '{\it{R}} (MPM)','Location','East')
xlabel('Incidence Angle (degrees)');
grid on


fig = figure(10); clf;
fig.Position = [100 100 1200 300];
subplot(121);
plot(thetas, T, thetas(1:3:end), mT(1:3:end),'x',thetas, R, thetas(1:3:end), mR(1:3:end),'o');
ylabel('{\it{R}} and {\it{T}}');
text(5, 0.5, ' TE_z ','EdgeColor', 'black')
legend('{\it{T}} (SIM)', '{\it{T}} (MPM)','{\it{R}} (SIM)', '{\it{R}} (MPM)','Location','NorthEast')
xlabel('Incidence Angle (degrees)');
text(-14, 1.05, '(a)', 'FontSize',24);
grid on
%
subplot(122);
plot(thetas, Ttm, thetas(1:3:end), mTtm(1:3:end),'x',thetas, Rtm, thetas(1:3:end), mRtm(1:3:end),'o');
ylabel('{\it{R}} and {\it{T}}');
text(5, 0.5, ' TM_z ','EdgeColor', 'black')
text(-14, 1.05, '(b)', 'FontSize',24);
legend('{\it{T}} (SIM)', '{\it{T}} (MPM)','{\it{R}} (SIM)', '{\it{R}} (MPM)','Location','East')
xlabel('Incidence Angle (degrees)');
grid on
