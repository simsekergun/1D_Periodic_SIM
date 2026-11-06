%% Accuracy vs sampling density N
%  All lengths are in nm.  
clear; 
set(0, 'DefaultFigureCreateFcn', @(h,~) set(h, 'Renderer', 'painters'))

%% Convergence vs N for 5 gratings, TE and TM
cases = { 'C1',   200, 633, eps_bk7(633),  500, 15;
          'C2',  300, 633, eps_bk7(633), 1500, 20;
          'C3',    100, 633, eps_gold(633), 500, 15;
          'C4',   150, 633, eps_gold(633),1200, 30;
          'C5',     100, 633, eps_gold(633), 220,  0};

Ns = [8 16 24 32 48 64 96 128 192 256];

pols = {'TE','TM'};
labelss = {'(a)', '(b)'};
err = nan(size(cases,1), numel(Ns), 2);  
res = err;
for c = 1:size(cases,1)
    [nm,a,lam,ep,d,thd] = cases{c,:};
    for ip = 1:2
        o2 = struct('pol',pols{ip},'thetai',thd*pi/180);
        ref = mpm1d(a, lam, ep, d, setfield(o2,'Lmax',60));
        for k = 1:numel(Ns)
            o = sim1dperiodic(struct('type','circle','a',a), lam, ep, Ns(k), d, o2);
            err(c,k,ip) = max([abs(o.Rn - ref.Rn), abs(o.Tn - ref.Tn)]);
            res(c,k,ip) = abs(o.energy_residual);
        end
        fprintf('%-28s %s R=%.6f T=%.6f A=%.6f | ', nm, pols{ip}, ref.R, ref.T, ref.A);
        fprintf('%8.0e', err(c,:,ip)); fprintf('\n');
    end
end

fig = figure(10); clf;
fig.Position = [100 100 1300  500];
for ip = 1:2
    subplot(1,2,ip); 
    semilogy(Ns, squeeze(err(:,:,ip)).', 'o-');     
    grid on;
    xlabel('{\it{N}}'); 
    ylabel('max |\Delta{\it{R}}_n|,|\Delta{\it{T}}_n|'); 
    text(110, 1e-2, [pols{ip} '_z'],'EdgeColor', 'black');
    xlim([0 260]);
    text(-55, 2, labelss{ip},'FontSize',24)
    legend(cases(:,1), 'Location','northeast');
end

fig = figure(11); clf;
fig.Position = [100 100 700  800];
for ip = 1:2
    subplot(2,1,ip); 
    semilogy(Ns, squeeze(err(:,:,ip)).', 'o-'); 
    grid on;
    xlabel('{\it{N}}'); 
    ylabel('max |\Delta{\it{R}}_n|,|\Delta{\it{T}}_n|'); 
    text(-40, 2, labelss{ip},'FontSize',24)
    text(110, 1e-2, [pols{ip} '_z'],'EdgeColor', 'black');
    xlim([0 260]);
    legend(cases(:,1), 'Location','northeast');
end
