function geom = rounded_rect_geom(A, B, p, M)
% Analytic rounded rectangle, equal-arclength parametrized,
%                   for sim1dperiodic / sim2d (geom.type = 'custom').
%
%   geom = rounded_rect_geom(A, B, p)       A, B: half-width/half-height
%   geom = rounded_rect_geom(A, B, p, M)    M: FFT grid for the arclength map
%
%   Boundary: polar superellipse rho(t) = ((cos t/A)^p + (sin t/B)^p)^(-1/p),
%   p even (analytic), evaluated in log form so that large p and large
%   aspect ratios do not overflow.  The curve is then reparametrized by
%   arclength, sigma = 2*pi*s/L, so that the N equispaced SIM nodes are
%   equispaced along the boundary (J = L/(2*pi)), which is important for
%   thin bars and sharp corners.  r, r', r'' are returned analytically
%   (chain rule through t(sigma)); t(sigma) is found by Newton iteration on
%   s(t) = sigma*L/(2*pi), with s(t) from the Fourier series of J(t).
%
%   Corner radius of curvature decreases with p (e.g. A = 1500 nm:
%   B = 300 nm, p = 46 gives about 15 nm).
%
% Ergun Simsek, UMBC

if nargin < 4, M = 2^14; end
if mod(p,2) ~= 0, error('rounded_rect_geom:p','p must be even'); end
t  = 2*pi*(0:M-1)/M;
[~, dr] = se_curve(t, A, B, p);
J  = hypot(dr(1,:), dr(2,:));
Jh = fft(J)/M;
if abs(Jh(M/2+1)) > 1e-12*abs(Jh(1))
    warning('rounded_rect_geom:M','Fourier tail of J not resolved; increase M');
end
L  = 2*pi*real(Jh(1));
k  = [0:M/2-1, -M/2:-1];
nz = k ~= 0;
ck = zeros(1,M);  
ck(nz) = Jh(nz)./(1i*k(nz));
sg = real(Jh(1))*t + real(ifft(ck)*M) - real(sum(ck));   % s(t) on the grid
sg(1) = 0;                                % exact start (FFT round-off)
D = struct('A',A,'B',B,'p',p,'L',L,'J0',real(Jh(1)), ...
           'sgrid',[sg L],'tgrid',[t 2*pi],'kk',k(nz),'cc',ck(nz).');
geom = struct('type','custom','A',A,'B',B,'p',p,'L',L);
geom.fun = @(sig) al_eval(sig, D);      % plain handle + data: safe in parfor
end

function [r, drs, ddrs] = al_eval(sig, D)
% equal-arclength evaluation: t(sigma) by Newton on s(t) = sigma*L/(2*pi)
sig = sig(:).';
tgt = mod(sig, 2*pi)*D.L/(2*pi);
tq  = interp1(D.sgrid, D.tgrid, tgt, 'linear', 'extrap');  % s -> t start
tq  = min(max(tq, 0), 2*pi);
for it = 1:6
    [~, d1] = se_curve(tq, D.A, D.B, D.p);
    Jq = hypot(d1(1,:), d1(2,:));
    sq = zeros(size(tq));
    for i0 = 1:256:numel(tq)                        % chunked Fourier series
        ii = i0:min(i0+255, numel(tq));
        sq(ii) = D.J0*tq(ii) + real((exp(1i*tq(ii).'*D.kk) - 1)*D.cc).';
    end
    f  = sq - tgt;
    tq = tq - f./Jq;
    if max(abs(f)) < 1e-12*D.L, break; end
end
if any(~isfinite(tq))
    error('rounded_rect_geom:nan','arclength inversion failed');
end
[r, d1, d2] = se_curve(tq, D.A, D.B, D.p);
Jq  = hypot(d1(1,:), d1(2,:));
tp  = D.L/(2*pi)./Jq;                                % dt/dsigma
Jp  = (d1(1,:).*d2(1,:) + d1(2,:).*d2(2,:))./Jq;     % dJ/dt
tpp = -tp.^2.*Jp./Jq;                                % d2t/dsigma2
drs  = d1.*tp;
ddrs = d2.*tp.^2 + d1.*tpp;
end

function [r, dr, ddr] = se_curve(t, A, B, p)
% polar superellipse and its t-derivatives, overflow-safe
c = cos(t);  s = sin(t);
l1 = p*log(abs(c)/A);  l2 = p*log(abs(s)/B);
lm = max(l1, l2);
logq = lm + log(exp(l1 - lm) + exp(l2 - lm));
T1 = exp(l1 - logq);  T2 = exp(l2 - logq);     % weights, T1 + T2 = 1
a1 = zeros(size(t));  
a2 = a1;  
b1 = a1;  
b2 = a1;
m1 = T1 > 0;  
m2 = T2 > 0;
a1(m1) = T1(m1).*s(m1)./c(m1);     
b1(m1) = T1(m1).*(s(m1)./c(m1)).^2;
a2(m2) = T2(m2).*c(m2)./s(m2);     
b2(m2) = T2(m2).*(c(m2)./s(m2)).^2;
Q1 = p*(-a1 + a2);                  % q'/q
Q2 = p*((p-1)*(b1 + b2) - 1);       % q''/q
rho  = exp(-logq/p);
rho1 = rho.*(-Q1/p);
rho2 = rho.*((1/p)*(1/p+1)*Q1.^2 - Q2/p);
r   = [rho.*c;  rho.*s];
dr  = [rho1.*c - rho.*s;  rho1.*s + rho.*c];
ddr = [(rho2 - rho).*c - 2*rho1.*s;  (rho2 - rho).*s + 2*rho1.*c];
end
