function [Sig, nvec, info] = lattice_sums_1d(k, beta, d, nmax, rho0f)
%LATTICE_SUMS_1D  Lattice sums of a 1-D array (period d along x)
%   Sigma_n = sum_{p~=0} exp(i beta p d) H_n^(1)(k|p|d) exp(i n arg(-p d x^)),
%           = sum_{p>=1} H_n^(1)(k p d) [ (-1)^n e^{i beta p d} + e^{-i beta p d} ],
%   n = -nmax..nmax, obtained from the Fourier coefficients of the Ewald
%   regular part on circles of radius rho = rho0f*d:
%   Gtilde(rho,phi) = (i/4) sum_m Sigma_{-m} J_m(k rho) e^{i m phi},  rho < d.
%
%   [Sig, nvec, info] = lattice_sums_1d(k, beta, d, nmax)
%   [Sig, nvec, info] = lattice_sums_1d(k, beta, d, nmax, rho0f)
%
%   rho0f  vector of radius fractions, 0 < rho0f < 1 (default [0.95 0.85 0.75]).
%          For each order m the coefficient is taken from the circle on which
%          |J_m(k rho)| is largest. A single radius fails where J_m(k rho) = 0
%          (e.g. d = 400 nm, rho = 0.95 d: J_1 = 0 at 623.12 nm, J_2 = 0 at
%          464.91 nm), because the Fourier coefficient is divided by J_m.
%          Passing a scalar rho0f reproduces the original single-radius method.
%   info   .rho0f, .radius_index (circle used for each nvec), .absJ (|J_m| used)
if nargin < 5 || isempty(rho0f), rho0f = [0.95 0.85 0.75]; end
if any(rho0f <= 0 | rho0f >= 1)
    error('lattice_sums_1d:radius', 'rho0f must lie in (0,1): the expansion holds for rho < d');
end
M  = 2^ceil(log2(max(2048, 4*nmax + 800)));
ph = (0:M-1)*2*pi/M;
m  = -nmax:nmax;
Sm = zeros(size(m));  
Jbest = zeros(size(m));  
ibest = zeros(size(m));
for ir = 1:numel(rho0f)
    rho0 = rho0f(ir)*d;
    Gt = greenper1d(rho0*cos(ph), rho0*sin(ph), k, beta, d);
    c  = fft(Gt)/M;                     % coefficient of e^{i m phi}
    J  = besselj(m, k*rho0);
    use = abs(J) > Jbest;               % keep the circle with the largest |J_m|
    Sm(use) = (4/1i)*c(mod(m(use),M)+1)./J(use);   % = Sigma_{-m}
    Jbest(use) = abs(J(use));  
    ibest(use) = ir;
end
nvec = -m;
[nvec, ix] = sort(nvec);  
Sig = Sm(ix);
info = struct('rho0f', rho0f, 'radius_index', ibest(ix), 'absJ', Jbest(ix));
end
