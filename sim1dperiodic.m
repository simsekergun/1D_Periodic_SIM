function out = sim1dperiodic(geom, lam, epsr, N, d, opts)
% Spectral integral method (SIM) for a 1-D periodic array
%   (grating) of identical homogeneous cylinders: reflection, transmission
%   and absorption of a plane wave, TE_z or TM_z, complex permittivity.
%
%   out = sim1dperiodic(geom, lam, epsr, N, d)
%   out = sim1dperiodic(geom, lam, epsr, N, d, opts)
%
%   Array: cylinders centred at (p d, 0), p integer, period d along x.
%   Incident unit plane wave from above (y > 0), angle thetai from the -y
%   axis (positive towards +x):
%       U_inc = exp(i(beta x - gamma0 y)),  beta = kb sin(thetai),
%                                           gamma0 = kb cos(thetai).
%   exp(-i w t) convention, Im(epsr) > 0 for lossy media.
%
%   geom, lam, epsr, N : as in sim2d (geom.type circle/ellipse/pert/rose/custom)
%   d     period (same unit as lam); the object must fit in |x| < d/2
%   opts  .pol ('TE'|'TM'), .thetai (rad, 0), .epsb (1), .Lcut (8), .Lon (18)
%         .H (1.5) Ewald splitting parameter (see greenper1d)
%
%   OUTPUT (in addition to the boundary data of sim2d: theta, r, Ja, nrm, U, dUn ...)
%     orders  propagating Floquet orders n (beta_n = beta + 2 pi n/d, |beta_n| < kb)
%     Rn, Tn  complex amplitudes: above  U_s = sum Rn exp(i(beta_n x + gamma_n y))
%                                 below  U   = sum Tn exp(i(beta_n x - gamma_n y))
%     R, T, A reflectance, transmittance, absorptance (R + T + A = 1)
%     energy_residual = 1 - R - T - A
%
%   Exterior kernel: G_p = (i/4)H0(kb R) + Gtilde, Gtilde from greenper1d
%   (smooth on the reference cylinder because 2a < d), so the SIM singularity
%   subtraction of sim2d is unchanged.  Interior kernel: free-space G_1.
%   Boundary conditions: TE dU/dn|in = (epsr/epsb) dU/dn|out; TM: equal.

if nargin < 6, opts = struct(); end
pol    = upper(getopt(opts,'pol','TE'));
thetai = getopt(opts,'thetai',0);
epsb   = getopt(opts,'epsb',1);
Lcut   = getopt(opts,'Lcut',8);
Lon    = getopt(opts,'Lon',18);
Hew    = getopt(opts,'H',1.5);
if mod(N,2)~=0, error('sim1dperiodic:N','N must be even'); end

% wavenumbers
k0 = 2*pi/lam;  
kb = k0*sqrt(epsb);  
k1 = k0*sqrt(epsr);
switch pol
    case 'TE', rho = epsr/epsb;
    case 'TM', rho = 1;
    otherwise, error('sim1dperiodic:pol','pol must be TE or TM');
end
beta = kb*sin(thetai);  
gam0 = kb*cos(thetai);

ii = 1:N;  
th = (ii-1)*2*pi/N;
[rr,dr,ddr] = curve(geom, th);
Ja  = sqrt(dr(1,:).^2 + dr(2,:).^2);
kap = (dr(1,:).*ddr(2,:) - dr(2,:).*ddr(1,:))./Ja.^3;
nrm = [dr(2,:); -dr(1,:)]./[Ja; Ja];
if 2*max(abs(rr(1,:))) >= d
    error('sim1dperiodic:fit','object does not fit in one period (need max|x| < d/2)');
end

tw = [Inf, Inf];
Rmax = max(max(hypot(rr(1,:).' - rr(1,:), rr(2,:).' - rr(2,:))));
if isfinite(Lon) && imag(k1)*Rmax > Lon
    tw(2) = Lcut/(imag(k1)*max(Ja));
end

%---- regular part of the periodic Green's function, all node pairs --------
Xm = rr(1,:).' - rr(1,:);           % X(j,i) = x_test(j) - x_src(i)
Ym = rr(2,:).' - rr(2,:);
[Gt, Gtx, Gty, ginfo] = greenper1d(Xm, Ym, kb, beta, d, struct('H',Hew));
Gtn = -(Gtx.*nrm(1,:) + Gty.*nrm(2,:));     % d Gtilde / d n'  (source normal)

%---- incident field --------------------------------------------------------
inc = exp(1i*(beta*rr(1,:) - gam0*rr(2,:)));
ff  = (fft(inc, N)/N).';              % e^{+ikt} convention, consistent with matrix rows

mf = 2*N;  
nn = 1:N-1;  
a0 = -1.0./(2.0*ii);
A = zeros(N); 
B = A; 
C = A; 
D = A; 
EE = A;
kk = [kb, k1];
for j = 1:N
    tm = th(j);  etm = exp(1i*(ii-1)*tm);  rt = rr(:,j);
    aj = zeros(1,mf);
    aj(2:N)    = a0(nn).*exp(1i*nn*tm);
    aj(N+1:mf) = a0(N:-1:1).*exp(-1i*(N:-1:1)*tm);
    ajk = fft(aj, mf);
    for s = 1:2
        k = kk(s);
        [g,h,p0,q0] = kernel_row(k, th, rr, Ja, nrm, kap, tm, rt, tw(s), j);
        if s == 1                                  % periodic exterior kernel
            g = g + ifft(Gt(j,:).*Ja, N);
            h = h + ifft(Gtn(j,:).*Ja, N);
        end
        u = padconv(p0, ajk, N);  
        v = padconv(q0, ajk, N);
        hk = -2*(2*pi*h + k.*v);   gk = 2*(-u + 2*pi*g);
        if s == 1, A(j,:) = hk + etm;  B(j,:) = gk;
        else,      C(j,:) = hk - etm;  D(j,:) = gk;  end
    end
    EE(j,:) = etm;
end
X   = (rho*D)\C;
phi = (A - B*X)\(EE*(2*ff));
psi = -X*phi;
U   = N*ifft(phi.', N);   
dUn = N*ifft(psi.', N);    % U(t_m) = sum_k phi_k e^{+ik t_m}

out = struct('theta',th,'r',rr,'Ja',Ja,'nrm',nrm,'U',U,'dUn',dUn, ...
    'dUn_in',rho*dUn,'Uinc',inc,'pol',pol,'lam',lam,'epsr',epsr,'epsb',epsb, ...
    'k0',k0,'kb',kb,'k1',k1,'beta',beta,'gamma0',gam0,'d',d,'thetai',thetai, ...
    'N',N,'tw',tw,'ewald',ginfo,'geom',geom);
out = floquet_rt(out);
end

%=========================================================================
function o = floquet_rt(o)
% Floquet reflection/transmission amplitudes from the boundary data
kb = o.kb;  
beta = o.beta;  
d = o.d;  
w = (2*pi/o.N)*o.Ja;
x = o.r(1,:);  
y = o.r(2,:);  
nx = o.nrm(1,:);  
ny = o.nrm(2,:);
nlo = floor((-kb-beta)*d/(2*pi));  
nhi = ceil((kb-beta)*d/(2*pi));
ords = nlo:nhi;  
bn = beta + 2*pi*ords/d;  
keep = abs(bn) < kb;
ords = ords(keep);  
bn = bn(keep);  
gn = sqrt(kb^2 - bn.^2);
Rn = zeros(size(ords));  
Tn = Rn;
for m = 1:numel(ords)
    Rn(m) = (1i/(2*d*gn(m)))*sum((-1i*(bn(m)*nx + gn(m)*ny).*o.U - o.dUn) ...
            .*exp(-1i*(bn(m)*x + gn(m)*y)).*w);
    Tn(m) = (ords(m)==0) + (1i/(2*d*gn(m)))*sum((-1i*(bn(m)*nx - gn(m)*ny).*o.U - o.dUn) ...
            .*exp(-1i*(bn(m)*x - gn(m)*y)).*w);
end
o.orders = ords;  
o.beta_n = bn;  
o.gamma_n = gn;  
o.Rn = Rn;  
o.Tn = Tn;
o.R = sum(gn/o.gamma0.*abs(Rn).^2);
o.T = sum(gn/o.gamma0.*abs(Tn).^2);
Cabs = -(1/kb)*sum(imag(conj(o.U).*o.dUn).*w);
o.A = Cabs/(d*cos(o.thetai));
o.energy_residual = 1 - o.R - o.T - o.A;
end

%=========================================================================
% The following helpers are identical to those in sim2d.m
function [g,h,p,q] = kernel_row(k, th, rr, Ja, nrm, kap, tm, rt, tw, jself)
euler = 0.57721566490153288;
d    = th - tm;
self = false(size(th)); 
self(jself) = true;
Rx = rr(1,:) - rt(1);  
Ry = rr(2,:) - rt(2);
R  = sqrt(Rx.^2 + Ry.^2);  R(self) = 1;
RnR = (Rx.*nrm(1,:) + Ry.*nrm(2,:))./R;
x  = k*R;
xl = zeros(size(th));  
xl(~self) = log(2*abs(sin(d(~self)/2)));
chi = cutoff(d, tw);  
act = chi > 0;
J0 = zeros(size(th));  
J1 = J0;
J0(act) = besselj(0, x(act));  
J1(act) = besselj(1, x(act));
H0 = besselh(0, 1, x);  
H1 = besselh(1, 1, x);
p = J0.*chi;          
p(self) = 1;
q = J1.*RnR.*chi;     
q(self) = 0;
g = (1i/4)*H0 + (1/(2*pi))*xl.*p;
g(self) = -(euler + log(k*Ja(self)/2))/(2*pi) + 1i/4;
h = -(1i/4)*k*RnR.*H1 - (k/(2*pi))*xl.*q;
h(self) = -kap(self)/(4*pi);
N = numel(th);
g = ifft(g.*Ja, N);  h = ifft(h.*Ja, N);
p = ifft(p.*Ja, N);  q = ifft(q.*Ja, N);
end

function u0 = padconv(p0, ajk, N)
mf = 2*N;  
p = zeros(1,mf);
p(1:N/2) = p0(1:N/2);  
p(mf-N/2+1:mf) = p0(N/2+1:N);
u = ifft(fft(p,mf).*ajk, mf);
u0 = [u(1:N/2), u(mf-N/2+1:mf)];
end

function c = cutoff(t, tw)
if ~isfinite(tw), c = ones(size(t)); return; end
t  = abs(angle(exp(1i*t)));  
s = min(t/tw, 1);
f1 = zeros(size(s));  
f2 = f1;
m1 = s > 0;  f1(m1) = exp(-1./s(m1));
m2 = s < 1;  f2(m2) = exp(-1./(1-s(m2)));
c  = f2./(f1 + f2);
end

function [r,dr,ddr] = curve(geom, t)
switch lower(geom.type)
    case 'circle'
        a = geom.a;
        r = [a*cos(t); a*sin(t)];  
        dr = [-a*sin(t); a*cos(t)];  
        ddr = -r;
    case 'ellipse'
        a = geom.a; b = geom.b;
        r = [a*cos(t); b*sin(t)];  
        dr = [-a*sin(t); b*cos(t)];  
        ddr = -r;
    case 'pert'
        a = geom.a; b = geom.b;  s = t + 0.05*sin(2*t);  f = 1 + 0.1*cos(2*t);
        r   = [ a*cos(s);  b*sin(s)];
        dr  = [-a*sin(s).*f;  b*cos(s).*f];
        ddr = [-r(1,:).*f.^2 + 0.2*a*sin(s).*sin(2*t); -r(2,:).*f.^2 - 0.2*b*cos(s).*sin(2*t)];
    case 'rose'
        a = geom.a; 
        dd = geom.d; 
        m = geom.m;
        rho = a*(1 + dd*cos(m*t));  
        drho = -a*dd*m*sin(m*t);  
        d2 = -a*dd*m^2*cos(m*t);
        r   = [rho.*cos(t); rho.*sin(t)];
        dr  = [drho.*cos(t) - rho.*sin(t); drho.*sin(t) + rho.*cos(t)];
        ddr = [(d2 - rho).*cos(t) - 2*drho.*sin(t); (d2 - rho).*sin(t) + 2*drho.*cos(t)];
    case 'custom'
        [r,dr,ddr] = geom.fun(t);
    otherwise
        error('sim1dperiodic:geom','unknown geometry type %s', geom.type);
end
end

function v = getopt(s, f, d)
if isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
