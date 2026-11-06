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

a = 200;        % radius
geom = struct('type','circle','a',a);    
d    = 800;                                         % period (nm)
lams = 600:900;         % wavelengths
theta_i = 10*pi/180;        % incidence angles (radian)
Nsim = 128;

R = zeros(length(lams),1);
T = R;
A = R;
res = R;
mR = R;
mT = R;
mA = R;
Rtm = R;
Ttm = R;
Atm = R;
restm = R;
mRtm = R;
mTtm = R;
mAtm = R;
kt = 1;
parfor k = 1:length(lams)
    ep = eps_bk7(lams(k));
    m = mpm1d(a, lams(k), ep, d, struct('pol','TE','thetai',theta_i));
    o = sim1dperiodic(geom, lams(k), ep, Nsim, d, struct('pol','TE','thetai',theta_i));
   mTM = mpm1d(a, lams(k), ep, d, struct('pol','TM','thetai',theta_i));
    oTM = sim1dperiodic(geom, lams(k), ep, Nsim, d, struct('pol','TM','thetai',theta_i));
  
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

fig = figure(9); clf;
fig.Position = [100 100 500 700];
subplot(211);
plot(lams, T, lams(1:10:end), mT(1:10:end),'x',lams, R, lams(1:10:end), mR(1:10:end),'o');
ylabel('{\it{R}} and {\it{T}}');
text(610, 0.5, ' TE_z ','EdgeColor', 'black')
legend('{\it{T}} (SIM)', '{\it{T}} (MPM)','{\it{R}} (SIM)', '{\it{R}} (MPM)','Location','East')
xlabel('Wavelength (nm)');
text(555, 1.05, '(a)', 'FontSize',24);
grid on
subplot(212);
plot(lams, Ttm, lams(1:10:end), mTtm(1:10:end),'x',lams, Rtm, lams(1:10:end), mRtm(1:10:end),'o');
ylabel('{\it{R}} and {\it{T}}');
text(610, 0.5, ' TM_z ','EdgeColor', 'black')
text(555, 1.05, '(b)', 'FontSize',24);
legend('{\it{T}} (SIM)', '{\it{T}} (MPM)','{\it{R}} (SIM)', '{\it{R}} (MPM)','Location','East')
xlabel('Wavelength (nm)');
grid on


fig = figure(10); clf;
fig.Position = [100 100 1200 300];
subplot(121);
plot(lams, T, lams(1:10:end), mT(1:10:end),'x',lams, R, lams(1:10:end), mR(1:10:end),'o');
ylabel('{\it{R}} and {\it{T}}');
text(610, 0.5, ' TE_z ','EdgeColor', 'black')
legend('{\it{T}} (SIM)', '{\it{T}} (MPM)','{\it{R}} (SIM)', '{\it{R}} (MPM)','Location','East')
xlabel('Wavelength (nm)');
text(555, 1.05, '(a)', 'FontSize',24);
grid on
subplot(122);
plot(lams, Ttm, lams(1:10:end), mTtm(1:10:end),'x',lams, Rtm, lams(1:10:end), mRtm(1:10:end),'o');
ylabel('{\it{R}} and {\it{T}}');
text(610, 0.5, ' TM_z ','EdgeColor', 'black')
text(555, 1.05, '(b)', 'FontSize',24);
legend('{\it{T}} (SIM)', '{\it{T}} (MPM)','{\it{R}} (SIM)', '{\it{R}} (MPM)','Location','East')
xlabel('Wavelength (nm)');
grid on
