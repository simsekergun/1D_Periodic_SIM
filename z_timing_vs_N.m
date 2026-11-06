% Compute time of sim1dperiodic versus the number of boundary
%   points N, for a gold and a BK7 cylinder grating.
%
%   Geometry: circular cylinders of radius a = lambda, gap lambda/2,
%   i.e. period d = 2a + lambda/2 = 2.5*lambda, in vacuum, normal incidence.
%   Each case is run once as a warm-up (not timed) and then nrep times;
%   the median and the min/max are reported and plotted.
%
% Ergun Simsek, UMBC

clear; clc;
set(0, 'DefaultFigureCreateFcn', @(h,~) set(h, 'Renderer', 'painters'))

lam  = 633;                      % wavelength (nm)
a    = lam;                        % radius (nm)
gap  = lam/2;                    % gap between neighbouring cylinders (nm)
d    = 2*a + gap;                % period (nm)
Ns   = [128 256 512 1024];
nrep = 10;
pol  = 'TE';                     % 'TE' or 'TM'
opts = struct('pol', pol, 'thetai', 0);
geom = struct('type', 'circle', 'a', a);

mats  = {'Gold', 'BK7'};
epss  = [eps_gold(lam), eps_bk7(lam)];

tAll  = nan(numel(mats), numel(Ns), nrep);
RTA   = nan(numel(mats), numel(Ns), 4);      % R, T, A, energy residual
fprintf('lambda = %g nm, a = %g nm, d = %g nm, %s, %d repetitions\n', ...
        lam, a, d, pol, nrep);

for im = 1:numel(mats)
    for k = 1:numel(Ns)
        N   = Ns(k);
        out = sim1dperiodic(geom, lam, epss(im), N, d, opts);   % warm-up
        for r = 1:nrep
            t0  = tic;
            out = sim1dperiodic(geom, lam, epss(im), N, d, opts);
            tAll(im,k,r) = toc(t0);
        end
        RTA(im,k,:) = [out.R, out.T, out.A, out.energy_residual];
        fprintf('%-5s N = %5d: median %8.4f s  (min %8.4f, max %8.4f)  R = %.10f  T = %.10f  res = %.1e\n', ...
                mats{im}, N, median(tAll(im,k,:)), min(tAll(im,k,:)), ...
                max(tAll(im,k,:)), out.R, out.T, out.energy_residual);
    end
end

tMed = median(tAll, 3);
tMin = min(tAll, [], 3);
tMax = max(tAll, [], 3);

% local and overall log-log slopes
for im = 1:numel(mats)
    p = polyfit(log(Ns), log(tMed(im,:)), 1);
    sl = diff(log(tMed(im,:))) ./ diff(log(Ns));
    fprintf('%-5s overall slope %.2f, local slopes %s\n', mats{im}, p(1), ...
            sprintf('%.2f ', sl));
end

save('timing_vs_N.mat', 'lam', 'a', 'd', 'Ns', 'nrep', 'pol', 'mats', ...
     'tAll', 'tMed', 'tMin', 'tMax', 'RTA');

% ------------------------------------------------------------------ plot
fig = figure('Color', 'w', 'Position', [200 200 520 300]);
mk  = {'o-', 's-'};
for im = 1:numel(mats)
    errorbar(Ns, tMed(im,:), tMed(im,:) - tMin(im,:), tMax(im,:) - tMed(im,:), ...
             mk{im}, 'LineWidth', 1.2, 'MarkerSize', 6); hold on;
end
% reference slopes anchored at the gold N = 512 point
Nr = [Ns(1) Ns(end)];
t0 = tMed(1, Ns == 512);
loglog(Nr, t0*(Nr/512).^1.5, 'k--', 'LineWidth', 1.2);
set(gca, 'XScale', 'log', 'YScale', 'log');
xticks(Ns); xlim([Ns(1)/1.2, Ns(end)*1.2]);
grid on; box on;
xlabel('{\itN}');
ylabel('Compute time (s)');
legend([mats, {'\propto {\itN}^{1.5}'}], 'Location', 'NorthWest');
