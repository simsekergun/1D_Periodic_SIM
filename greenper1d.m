function [G, Gx, Gy, info] = greenper1d(X, Y, k, beta, d, opts)
%GREENPER1D  1-D quasi-periodic Green's function of the 2-D Helmholtz equation
%            (Ewald representation) and its gradient.
%
%   G_p(X,Y) = (i/4) sum_p exp(i beta p d) H0^(1)( k |(X - p d, Y)| ),
%   (lattice along x, period d, Bloch wavenumber beta, exp(-i w t)).
%
%   [G,Gx,Gy] = greenper1d(X,Y,k,beta,d)               regular part (default)
%       Gtilde = G_p - (i/4) H0^(1)(k R),  R = sqrt(X^2+Y^2); smooth for |X| < d,
%       including the limit at R = 0.
%   [...] = greenper1d(X,Y,k,beta,d,struct('regular',false))   full G_p
%   opts.H   (default 1.5) splitting E = max(sqrt(pi)/d, k/(2H))
%   opts.tol (default 6.5) truncation: erfc / exp(-z) arguments > tol^2
%
%   Gx, Gy are derivatives with respect to X and Y (X = x - x').
%   Spectral part (n sums over Floquet orders, beta_n = beta + 2 pi n/d,
%   kappa_n = sqrt(beta_n^2 - k^2), Re >= 0, = -i gamma_n for propagating):
%     G_spec = 1/(4d) sum_n e^{i beta_n X}/kappa_n [P(kappa_n,Y) + P(kappa_n,-Y)]
%     P(kappa,Y) = e^{kappa Y} erfc(kappa/(2E) + Y E)
%   Spatial part:
%     G_spat = 1/(4 pi) sum_p e^{i beta p d} sum_q (k/2E)^{2q}/q! E_{q+1}(R_p^2 E^2)
%   Wood anomalies (kappa_n = 0) are singular and must be avoided.

if nargin < 6, opts = struct(); end
regular = getopt(opts,'regular',true);
H   = getopt(opts,'H',1.5);
tol = getopt(opts,'tol',6.5);
euler = 0.57721566490153288;

E = max(sqrt(pi)/d, k/(2*H));
ymax = max(abs(Y(:)));  
xmax = max(abs(X(:)));
G = zeros(size(X));  
Gx = G;  
Gy = G;

%---- spectral part ------------------------------------------------------
nmax = ceil((2*E*(tol + ymax*E) + k)*d/(2*pi)) + 2;
kapmin = Inf;
for n = -nmax:nmax
    bn  = beta + 2*pi*n/d;
    kap = -1i*sqrt(k^2 - bn^2);
    if real(kap) < 0, kap = -kap; end
    kapmin = min(kapmin, abs(kap));
    ph = exp(1i*bn*X);
    Pp = Pfun(kap,  Y, E);
    Pm = Pfun(kap, -Y, E);
    G  = G  + ph.*(Pp + Pm)/(4*d*kap);
    Gx = Gx + 1i*bn*ph.*(Pp + Pm)/(4*d*kap);
    Gy = Gy + ph.*(Pp - Pm)/(4*d);
end
if kapmin*d < 1e-6
    warning('greenper1d:wood','Close to a Rayleigh-Wood anomaly (|kappa_n| d = %.1e).', kapmin*d);
end

%---- spatial part -------------------------------------------------------
r2 = (k/(2*E))^2;  
c = 1;  
q = 0;
while true
    q = q + 1;  
    c(q+1) = c(q)*r2/q;
    if c(q+1)/q < 1e-18 && q > 3, break; end
end
qmax = q;
pmax = ceil((tol/E + xmax)/d) + 1;
for p = -pmax:pmax
    Xp = X - p*d;  
    Rp2 = Xp.^2 + Y.^2;  
    z = Rp2*E^2;
    ph = exp(1i*beta*p*d);
    if p == 0 && regular
        m  = z > 0;
        zz = z;  
        zz(~m) = 1;
        [Sq, Sq1] = qsums(zz, c, qmax);
        R  = sqrt(zz)/E;   
        kR = k*R;
        val = Sq/(4*pi) - (1i/4)*besselh(0,1,kR);
        lim = log(k/(2*E))/(2*pi) + euler/(4*pi)+ sum(c(2:end)./(1:qmax))/(4*pi) - 1i/4;
        val(~m) = lim;
        fac = -(E^2/(2*pi))*Sq1 + (1i*k/4)*besselh(1,1,kR)./R;
        fac(~m) = 0;                          % smooth & radially symmetric
        G = G + val;  
        Gx = Gx + fac.*Xp;  
        Gy = Gy + fac.*Y;
    else
        if any(z(:) == 0), error('greenper1d:sing','field point on a source'); end
        [Sq, Sq1] = qsums(z, c, qmax);
        G  = G + ph*Sq/(4*pi);
        fac = -ph*(E^2/(2*pi))*Sq1;
        Gx = Gx + fac.*Xp;  
        Gy = Gy + fac.*Y;
    end
end
info = struct('E',E,'nmax',nmax,'pmax',pmax,'qmax',qmax,'kapmin',kapmin);
end

%-------------------------------------------------------------------------
function P = Pfun(kap, y, E)
% e^{kap y} erfc(kap/(2E) + y E), evaluated without overflow via Faddeeva w
a = kap/(2*E) + y*E;  
z = 1i*a;
pref = exp(-kap^2/(4*E^2) - (y*E).^2);
neg = imag(z) < 0;
zz = z;  
zz(neg) = -z(neg);
w = faddeeva_w(zz);
P = pref.*w;
P(neg) = 2*exp(kap*y(neg)) - pref(neg).*w(neg);
end

function [S, S1] = qsums(z, c, qmax)
% S = sum_q c_q E_{q+1}(z),  S1 = sum_q c_q E_q(z),  E_0 = e^{-z}/z
ez = exp(-z);
Eq  = ez./z;          % E_0
Eq1 = expint(z);      % E_1
S = c(1)*Eq1;  
S1 = c(1)*Eq;
for q = 1:qmax
    Enew = (ez - z.*Eq1)/q;            % E_{q+1}
    Eq = Eq1;  Eq1 = Enew;
    S  = S  + c(q+1)*Eq1;
    S1 = S1 + c(q+1)*Eq;
end
end

function v = getopt(s, f, d)
if isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
