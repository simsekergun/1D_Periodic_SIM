function M = mpm1d(a, lam, epsr, d, opts)
%MPM1D  Independent reference: multipole (Rayleigh / lattice-sum) method for a
%       1-D periodic array of CIRCULAR cylinders (radius a, period d along x).
%
%   M = mpm1d(a, lam, epsr, d, opts)   opts: .pol, .thetai, .epsb, .Lmax
%   Same incidence and output conventions as sim1dperiodic (orders, Rn, Tn, R, T).
%
%   Scattered field  U_s = sum_l b_l sum_p e^{i beta p d} H_l(k rho_p) e^{i l phi_p}
%   Exciting field on cylinder 0: e = a_inc + L b,  L_{ml} = Sigma_{l-m}
%   Single-cylinder T-matrix: b_m = t_m e_m  (t_m = Mie coefficient a_m)
%   Solved in the scaled variable bh = H_l(ka) b_l for conditioning.
%   Floquet amplitudes:
%     R_n = 2/(d gamma_n) sum_l b_l (-i)^l ((beta_n + i gamma_n)/k)^l
%     T_n = delta_n0 + 2/(d gamma_n) sum_l b_l (-i)^l ((beta_n - i gamma_n)/k)^l
if nargin < 5, opts = struct(); end
pol    = upper(getopt(opts,'pol','TE'));
thetai = getopt(opts,'thetai',0);
epsb   = getopt(opts,'epsb',1);
k0 = 2*pi/lam;  
kb = k0*sqrt(epsb);
Lmax = getopt(opts,'Lmax', ceil(kb*a + 4*(kb*a)^(1/3) + 30));
beta = kb*sin(thetai);  
gam0 = kb*cos(thetai);

Mie = mie_cyl(a, lam, epsr, 0, struct('pol',pol,'epsb',epsb,'nmax',Lmax));
t = Mie.an;  
l = Mie.n;                          % l = -Lmax..Lmax
[Sig, nv] = lattice_sums_1d(kb, beta, d, 2*Lmax);
L = Sig(l.' - l + 2*Lmax + 1);                   % L(m,l) = Sigma_{l-m}
L = L.';                                         % rows m, columns l
Hl = besselh(l, 1, kb*a);
ainc = (1i.^l).*exp(-1i*l*(thetai - pi/2));
Amat = eye(numel(l)) - (t.*Hl).'.*L./Hl;
bh = Amat\((t.*Hl).*ainc).';
b  = bh.'./Hl;

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
    ep = ((bn(m) + 1i*gn(m))/kb).^l;  
    em = ((bn(m) - 1i*gn(m))/kb).^l;
    Rn(m) = 2/(d*gn(m))*sum(b.*(-1i).^l.*ep);
    Tn(m) = (ords(m)==0) + 2/(d*gn(m))*sum(b.*(-1i).^l.*em);
end
M = struct('orders',ords,'Rn',Rn,'Tn',Tn,'R',sum(gn/gam0.*abs(Rn).^2), ...
           'T',sum(gn/gam0.*abs(Tn).^2),'b',b,'l',l,'cond',cond(Amat));
M.A = 1 - M.R - M.T;
end

function v = getopt(s, f, d)
if isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
