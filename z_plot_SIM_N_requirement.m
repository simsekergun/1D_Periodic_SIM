clear;
set(0, 'DefaultFigureCreateFcn', @(h,~) set(h, 'Renderer', 'painters'))

a = 100;                % radius
lambda0 = 633;      % wavelength
Ns = [8 16 24 32 48 64 96 128 192 256 384];

ff = [0.4 0.6 0.8 0.9];   % filling fraction 2a/d
eg = nan(numel(ff), numel(Ns));


for f = 1:numel(ff)    
    d = 2*a/ff(f);      % periodicity 
    ref = mpm1d(a, lambda0, eps_gold(lambda0), d, struct('Lmax',60));
    for k = 1:numel(Ns)
        o = sim1dperiodic(struct('type','circle','a',a), lambda0, eps_gold(633), Ns(k), d);
        eg(f,k) = max([abs(o.Rn - ref.Rn), abs(o.Tn - ref.Tn)]);
    end
    fprintf('2a/d = %.2f: ', ff(f)); fprintf('%8.0e', eg(f,:)); fprintf('\n');
end


for iii = 1:4
    d = 2*a/ff(iii);
    g(iii) = d-2*a;
    ci = 2*pi*a;
    asd= find(eg(iii,:)<1e-10);
    Nmin(iii) = asd(1);
end


% What should be N for error < 1e-10 vs perimeter/gap 
tol   = 1e-10;
gap   = 2*a./ff - 2*a;            % gap between neighbouring cylinders (nm)
ratio = 2*pi*a ./ gap;            % perimeter / gap
Nreq  = nan(1, numel(ff));
for f = 1:numel(ff)
    idx = find(eg(f,:) < tol, 1, 'first');
    if ~isempty(idx) 
        Nreq(f) = Ns(idx); 
    end   % store N itself, not the index
end
disp(table(ff(:), gap(:), ratio(:), Nreq(:), 'VariableNames', {'ff','gap_nm','perim_over_gap','Nreq'}));

ok = ~isnan(Nreq);
p  = polyfit(ratio(ok), Nreq(ok), 1);
xr = linspace(0, max(ratio)*1.05, 100);


figure(16); clf; 
semilogy(Ns, eg,'-o'); 
grid on; 
legend('2a/d = 0.4','2a/d = 0.6','2a/d = 0.8','2a/d = 0.9'); 
xlabel('\it{N}'); 
ylabel('max |\Delta{\it{R}}_n|,|\Delta{\it{T}}_n|'); 

figure(17); clf;
plot(ratio(ok), Nreq(ok), 'o', xr, polyval(p, xr), '--'); grid on;
xlabel('Perimeter / gap');
ylabel('{\itN} required for \delta < 10^{-10}');
legend('SIM', sprintf('{\\itN} = %.1f\\times(2\\pi a/g) + %.0f', p(1), p(2)),'Location','NorthWest');


fig = figure(18); clf;
fig.Position = [100 100 600 800];
subplot(211); 
semilogy(Ns, eg,'-o'); 
grid on; 
legend('2a/d = 0.4','2a/d = 0.6','2a/d = 0.8','2a/d = 0.9'); 
xlabel('\it{N}'); 
ylabel('max |\Delta{\it{R}}_n|,|\Delta{\it{T}}_n|'); 

subplot(212); 
plot(ratio(ok), Nreq(ok), 'o', xr, polyval(p, xr), '--'); grid on;
xlabel('Perimeter / gap');
ylabel('{\itN} required for \delta < 10^{-10}');
legend('SIM', sprintf('{\\itN} = %.1f\\times(2\\pi a/g) + %.0f', p(1), p(2)),'Location','NorthWest');

fig = figure(19); clf;
fig.Position = [100 100 900 400];
subplot(121); 
semilogy(Ns, eg,'-o'); 
grid on; 
legend('2a/d = 0.4','2a/d = 0.6','2a/d = 0.8','2a/d = 0.9'); 
xlabel('\it{N}'); 
ylabel('max |\Delta{\it{R}}_n|,|\Delta{\it{T}}_n|'); 
text(-120, 2, '(a)','Fontsize',20)
text(430, 2, '(b)','Fontsize',20)
subplot(122); 
plot(ratio(ok), Nreq(ok), 'o', xr, polyval(p, xr), '--'); grid on;
xlabel('Perimeter / gap');
ylabel('{\itN} required for \delta < 10^{-10}');
legend('SIM', sprintf('{\\itN} = %.1f\\times(2\\pi a/g) + %.0f', p(1), p(2)),'Location','NorthWest');
