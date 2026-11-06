%% WOOD_BENCHMARK  
% SIM vs multipole near the (-1) Rayleigh-Wood anomaly
%   Gold circles a = 120 nm, d = 400 nm, theta_i = 30 deg -> lambda_RA = 600 nm.
clear;
set(0, 'DefaultFigureCreateFcn', @(h,~) set(h, 'Renderer', 'painters'))

a = 120;  % radius (nm)
d = 400;   % periodicity (nm)
th = 30*pi/180;     % incidence angle (radian)
N = 256;                % number of mesh points SIM 
dl = reshape([1;-1]*10.^(1:-1:-6), 1, []);
pols = {'TE','TM'};
fprintf('%4s %10s %20s %20s %10s\n','pol','dlam(nm)','R (SIM)','T (SIM)','|SIM-MPM|');
for ip = 1:2
    for k = 1:numel(dl)
        lam = 600 + dl(k); 
        ep = eps_gold(lam);
        o = sim1dperiodic(struct('type','circle','a',a), lam, ep, N, d, struct('pol',pols{ip},'thetai',th));
        m = mpm1d(a, lam, ep, d, struct('pol',pols{ip},'thetai',th));
        fprintf('%4s %10.0e %20.14f %20.14f %10.1e\n', pols{ip}, dl(k), o.R, o.T, max(abs([o.R-m.R, o.T-m.T])));
        if ip == 1
            lam_TE(k) = lam;
            R_TE(k) = o.R;
            T_TE(k) = o.T;
            error_TE(k) = max(abs([o.R-m.R, o.T-m.T]));
        else
            lam_TM(k) = lam;
            R_TM(k) = o.R;
            T_TM(k) = o.T;
            error_TM(k) = max(abs([o.R-m.R, o.T-m.T]));
        end
    end
end

figure(5); clf;
loglog(abs(lam_TE-600), error_TE,'o', abs(lam_TM-600), error_TM,'s','MarkerSize', 12);
xlabel('\delta\lambda (nm)');
ylabel('Difference between SIM and MPM');
legend('TE_z','TM_z')
