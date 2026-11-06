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
lams = 450:10:700;      % wavelengths (nm)
thetas = linspace(0, 75,26)*pi/180; % incidence angles (nm)
[x, y] = meshgrid(lams, thetas);        
Nsim = 128;     % Mesh sampling number

R = zeros(length(lams), length(thetas)); 
T = R; 
A = R; 
res = R;
for k = 1:numel(lams)
    parfor kt = 1:length(thetas)
        o = sim1dperiodic(geom, lams(k), eps_gold(lams(k)), Nsim, d, struct('pol','TM','thetai',thetas(kt)));
        R(k, kt) = o.R; 
        T(k, kt) = o.T; 
        A(k, kt) = o.A; 
        res(k, kt) = o.energy_residual;
    end
end

% save('results_TM.mat');

figure(6);
subplot(221); pcolor(x.', y.'*180/pi, R); 
colorbar; 
title('R'); 
shading interp;
xlabel('Wavelength (nm)');
ylabel('\theta_i (degrees)');
%
subplot(222); 
pcolor(x.', y.'*180/pi,T); 
colorbar; 
title('T'); 
shading interp;
xlabel('Wavelength (nm)');
ylabel('\theta_i (degrees)');
%
subplot(223); 
pcolor(x.', y.'*180/pi,A); 
colorbar;  
title('A'); 
shading interp;
xlabel('Wavelength (nm)');
ylabel('\theta_i (degrees)');
%
subplot(224); 
pcolor(x.', y.'*180/pi,res); 
colorbar;  
title('Residual'); 
shading interp;
xlabel('Wavelength (nm)');
ylabel('\theta_i (degrees)');
% print -dpng figure_RTA_lambda_vs_thetas_TM