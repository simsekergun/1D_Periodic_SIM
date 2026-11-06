%VALIDATE_SIM1DPERIODIC  Accuracy of sim1dperiodic vs sampling density N, with the
%   independent multipole (lattice-sum) solution mpm1d as reference.
%   Lengths in nm.  Expected values (from the NumPy port) are in the comments.
clear; close all;
set(0, 'DefaultFigureCreateFcn', @(h,~) set(h, 'Renderer', 'painters'))

%% 0. Green's function check: Ewald vs direct spectral sum at |Y| >= 0.15 d
d = 1; k = 2*pi/0.7; beta = k*sin(20*pi/180);
X = rand(1,200)*2-1;  Y = (rand(1,200)*0.6+0.15).*sign(rand(1,200)-0.5);
G = greenper1d(X, Y, k, beta, d, struct('regular',false));
Gs = zeros(size(X));
for n = -400:400
    bn = beta + 2*pi*n/d;  gn = sqrt(k^2 - bn^2);  if imag(gn) < 0, gn = -gn; end
    Gs = Gs + (1i/(2*d))*exp(1i*bn*X + 1i*gn*abs(Y))/gn;
end
fprintf('Ewald vs spectral sum: %.1e   (expected ~1e-15)\n', max(abs(G-Gs))/max(abs(Gs)));

%% 1. Convergence vs N for 7 gratings, TE and TM
cases = { 'C1',   200, 633, eps_bk7(633),  500, 15;
          'C2',  300, 633, eps_bk7(633), 1500, 20;
          'C3',    100, 633, eps_gold(633), 500, 15;
          'C4',   150, 633, eps_gold(633),1200, 30;
          'C5',     100, 633, eps_gold(633), 220,  0;
          'C6',   240, 633, eps_bk7(633),  500, 10;
          'C7',400, 800, eps_gold(800),1000,  5};
Ns = [8 16 24 32 48 64 96 128 192 256 384];
pols = {'TE','TM'};
err = nan(size(cases,1), numel(Ns), 2);  res = err;
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
% expected (TE): C1 7e-14 @64; C2 2e-15 @32; C3 9e-14 @24; C4 1e-13 @32;
%                C5 1e-11 @128, 2e-14 @192; C6 1e-13 @384; C7 8e-9 @384

Nnormalizer = 1;
% Nnormalizer = real(sqrt(vertcat(cases{:, 4}))).*2*pi.*vertcat(cases{:, 2})./vertcat(cases{:, 5});
NN = repmat(Ns, length(Nnormalizer), 1)./repmat(Nnormalizer,1,length(Ns));

labelss = {'(a)', '(b)'};
fig = figure(1); clf;
fig.Position = [100 100 1300  500];
for ip = 1:2
    subplot(1,2,ip); 
    semilogy(NN.', squeeze(err(:,:,ip)).', 'o-');     
    grid on;
    xlabel('{\it{N}}'); 
    ylabel('max |\Delta{\it{R}}_n|,|\Delta{\it{T}}_n|'); 
    text(210, 1e-2, [pols{ip} '_z'],'EdgeColor', 'black');
    %xlim([0 60]);
    text(-90, 2, labelss{ip},'FontSize',24)
    legend(cases(:,1), 'Location','northeast');
end
print -dpng ../docs/latex_source/sim1dperiodic/figure_SIM_accuracy_horiz

fig = figure(11); clf;
fig.Position = [100 100 700  800];
for ip = 1:2
    subplot(2,1,ip); 
    semilogy(NN.', squeeze(err(:,:,ip)).', 'o-'); 
    grid on;
    xlabel('{\it{N}}'); 
    ylabel('max |\Delta{\it{R}}_n|,|\Delta{\it{T}}_n|'); 
    text(-60, 2, labelss{ip},'FontSize',24)
    text(210, 1e-2, [pols{ip} '_z'],'EdgeColor', 'black');
    %xlim([0 60]);
    legend(cases(:,1), 'Location','northeast');
end
print -dpng ../docs/latex_source/sim1dperiodic/figure_SIM_accuracy_vertical
keyboard
%% 2. Gap study (gold a = 100 nm, 633 nm, TE, normal incidence)
fill = [0.4 0.6 0.8 0.9];  
eg = nan(numel(fill), numel(Ns));
for f = 1:numel(fill)
    a = 100; 
    d = 2*a/fill(f);
    ref = mpm1d(a, 633, eps_gold(633), d, struct('Lmax',60));
    for k = 1:numel(Ns)
        o = sim1dperiodic(struct('type','circle','a',a), 633, eps_gold(633), Ns(k), d);
        eg(f,k) = max([abs(o.Rn - ref.Rn), abs(o.Tn - ref.Tn)]);
    end
    fprintf('2a/d = %.2f: ', fill(f)); fprintf('%8.0e', eg(f,:)); fprintf('\n');
end
figure(16); clf; 
semilogy(Ns, eg,'-o'); 
grid on; 
legend('2a/d = 0.4','2a/d = 0.6','2a/d = 0.8','2a/d = 0.9'); 
xlabel('\it{N}'); 
ylabel('max |\Delta{\it{R}}_n|,|\Delta{\it{T}}_n|'); 
print -dpng ../docs/latex_source/sim1dperiodic/figure_SIM_gap_effect
% expected N for 1e-10: 24, 32, 64, 128  (rule N ~ 4.1*(2*pi*a/gap) + 14)

for iii = 1:4
    a = 100; 
    d = 2*a/fill(iii);
    g(iii) = d-2*a;
    ci = 2*pi*a;
    asd= find(eg(iii,:)<1e-10);
    Nmin(iii) = asd(1);
end

ok = ~isnan(Nreq);
p  = polyfit(ratio(ok), Nreq(ok), 1);
figure(17); clf;
plot(ratio(ok), Nreq(ok), 'o', ratio(ok), polyval(p, ratio(ok)), '--');
xlabel('Perimeter / gap');
ylabel('{\it N} required for \delta < 10^{-10}');
legend('SIM', '{\it{N}} = 4\times(2\pi a/g)+14', 'Location','NorthWest');
grid on;
print -dpng ../docs/latex_source/sim1dperiodic/figure_SIM_gap_N_requirement

%% 3. Spectrum: gold a=50, d=500, normal incidence, TE (N = 48) vs multipole
lams = 405:10:895;  Rs = zeros(size(lams)); Ts = Rs; As = Rs; Rm = Rs; Tm = Rs;
for k = 1:numel(lams)
    o = sim1dperiodic(struct('type','circle','a',50), lams(k), eps_gold(lams(k)), 48, 500);
    m = mpm1d(50, lams(k), eps_gold(lams(k)), 500, struct('Lmax',40));
    Rs(k) = o.R; Ts(k) = o.T; As(k) = o.A; Rm(k) = m.R; Tm(k) = m.T;
end
fprintf('gold grating spectrum: max|dR| = %.1e, max|dT| = %.1e  (expected ~5e-15)\n', ...
        max(abs(Rs-Rm)), max(abs(Ts-Tm)));
figure(2); clf; plot(lams, Rm, 'b-', lams, Rs, 'bo', lams, Tm, 'g-', lams, Ts, 'go', lams, As, 'rs');
xlabel('\lambda (nm)'); legend('R mpm','R SIM','T mpm','T SIM','A SIM');
