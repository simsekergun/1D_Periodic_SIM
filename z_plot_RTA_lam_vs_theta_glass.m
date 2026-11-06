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

geom = struct('type','circle','a',120);    % nm
d    = 400;                                         % period (nm)
lams = 600:5:900;                               % wavelengths (nm)
thetas = linspace(0, 75,16)*pi/180;     % incidence angles (radian)
[x, y] = meshgrid(lams, thetas);
Nsim = 128;                                     % SIM mesh number

R = zeros(length(lams), length(thetas)); 
T = R; 
A = R; 
res = R;
for k = 1:numel(lams)
    parfor kt = 1:length(thetas)
        o = sim1dperiodic(geom, lams(k), eps_bk7(lams(k)), Nsim, d, struct('pol','TE','thetai',thetas(kt)));
        R(k, kt) = o.R; 
        T(k, kt) = o.T; 
        A(k, kt) = o.A; 
        res(k, kt) = o.energy_residual;
    end
end


fig = figure(7);
fig.Position = [100 100 500 1200];
subplot(311); 
pcolor(x.', y.'*180/pi, R); 
colorbar; 
title('R'); 
shading interp;
xlabel('Wavelength (nm)');
ylabel('\theta_i (degrees)');
%
subplot(312); 
pcolor(x.', y.'*180/pi,T); 
colorbar; 
title('T'); 
shading interp;
xlabel('Wavelength (nm)');
ylabel('\theta_i (degrees)');
%
subplot(313); 
pcolor(x.', y.'*180/pi,A);
colorbar;  
title('A'); 
shading interp;
xlabel('Wavelength (nm)');
ylabel('\theta_i (degrees)');
% print -dpng figure_RTA_lambda_vs_thetas_glass