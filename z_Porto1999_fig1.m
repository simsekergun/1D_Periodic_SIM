%  Quick reproduction of Fig. 1 of 
%  J. A. Porto,F. J. Garcia-Vidal and J. B. Pendry, "Transmission resonances on 
%   metallic gratings with very narrow slits," PRL 83, 2845 (1999).
%
%   Free-standing gold lamellar grating in vacuum, d = 3.5 um, slit 0.5 um
%   (gold bars 3.0 um wide), heights h = 0.2 ... 4 um, normal incidence,
%   H parallel to the slits (TE_z here). Plotted: zero-order transmittance.
%
% Ergun Simsek, UMBC (November, 2026)
 
clear; 
clc;
d     = 3500;  % periodicity (nm)
aslit = 500;  % gap (nm)
A = (d - aslit)/2;   % width (nm)
hs    = [200 600 1200];           % grating heights (nm) 

% create a superellipse that is very similar to a rectangle
p     = 32;                                      % same exponent for all h
N     = 1024;                                    % boundary points
opts  = struct('pol','TE','thetai',0);

% coarse wavelength grid (nm), denser just above lambda = d
lam = [1000:150:3500, 3510:10:4200, 4300:100:5000, 5250:250:15000];
lam = lam(abs(lam - d) > 1 & abs(lam - d/2) > 1 & abs(lam - d/3) > 1);
epsAu = eps_gold_drude_ir(lam);

T0 = nan(numel(hs), numel(lam));
tAll = tic;
for ih = 1:numel(hs)
    geom = rounded_rect_geom(A, hs(ih)/2, p);
    fprintf('h = %4d nm: perimeter %.0f nm, N = %d, %d wavelengths ... ', ...
hs(ih), geom.L, N, numel(lam));
    tH = tic;
    Tk = nan(1, numel(lam));
    parfor il = 1:numel(lam)
        out    = sim1dperiodic(geom, lam(il), epsAu(il), N, d, opts);
        Tk(il) = abs(out.Tn(out.orders == 0))^2;
    end
    T0(ih,:) = Tk;
    fprintf('%.0f s\n', toc(tH));
    save('Porto1999_fig1_data.mat', 'lam', 'hs', 'p', 'N', 'T0');
end
fprintf('total %.1f min\n', toc(tAll)/60);

% plot results
fig = figure('Color','w','Position',[100 100 520 600]);
sty = {':','-','--'};
hold on;
for k = 1:length(hs)
    plot(lam/1000, T0(k,:), ['k' sty{k}], 'LineWidth', 1);
end
legend(compose('h = %.1f \\mum', hs/1000), 'Box','off', 'Location','northeast');
ylabel('Zero-order transmittance'); 
xlabel('Wavelength (\mum)');
grid on;
xlim([1 15]); 
ylim([0 1]); 
box on;