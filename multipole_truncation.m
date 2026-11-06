%MULTIPOLE_TRUNCATION  Which multipole orders are needed for the TE reflection band?
%   Gold circles a = 120 nm, d = 400 nm. 
%   For each theta, the R maximum on the
%   red side of the Rayleigh line, with the multipole method truncated at |n| <= Lmax.
%   Expected theta = 20 deg -> 
%   Lmax=1: peak pinned to lambda_RA, R = 0.475;
%   Lmax=2: +5.75 nm, R = 0.635; 
%   Lmax=3 and full: +7.674 nm, R = 0.619.
clear;
set(0, 'DefaultFigureCreateFcn', @(h,~) set(h, 'Renderer', 'painters'))

parpool_no = 8;

p = gcp('nocreate');
if isempty(p)
    parpool(parpool_no);
elseif p.NumWorkers ~= parpool_no
    delete(gcp('nocreate'));
    parpool(parpool_no);
end


a = 120; 
d = 400;
for thd = [14 20 30]
    th = thd*pi/180; 
    lra = d*(1 + sin(th)); 
    lam = lra + logspace(log10(0.003), log10(45), 300);
    for Lm = [1 2 3 4 0]
        R = zeros(size(lam));
        parfor k = 1:numel(lam)
            op = struct('pol','TE','thetai',th); 
            if Lm > 0, op.Lmax = Lm; end
            M = mpm1d(a, lam(k), eps_gold(lam(k)), d, op); 
            R(k) = M.R;
        end
        [Rp, i] = max(R);
        fprintf('theta %2d  Lmax %d : peak at lambda_RA + %.3f nm, R = %.3f\n', thd, Lm, lam(i)-lra, Rp);
    end
end
