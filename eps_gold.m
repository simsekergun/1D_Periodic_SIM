function e = eps_gold(lam_nm, table)
%EPS_GOLD  Relative permittivity of gold, exp(-i w t) convention (Im e > 0).
%   e = eps_gold(lam_nm)
%       Analytic Drude + two critical-point model fitted to Johnson & Christy:
%       P. G. Etchegoin, E. C. Le Ru, M. Meyer, J. Chem. Phys. 125, 164705
%       (2006); erratum 127, 189901 (2007).  Valid roughly 400-1000 nm.
%   e = eps_gold(lam_nm, table)
%       Interpolate user-supplied tabulated data instead,
%       table = [lambda_nm, n, k]  (e.g. Johnson & Christy from
%       refractiveindex.info);  e = (n + i k)^2.

lam_nm = double(lam_nm);
if nargin > 1 && ~isempty(table)
    n = interp1(table(:,1), table(:,2), lam_nm, 'pchip');
    k = interp1(table(:,1), table(:,3), lam_nm, 'pchip');
    e = (n + 1i*k).^2;
    return
end
einf = 1.53;  lp = 145;  gp = 17000;                  % nm
A  = [0.94, 1.36];   ph = [-pi/4, -pi/4];
li = [468, 331];     gi = [2300, 940];                % nm
l  = lam_nm;
e  = einf - 1./(lp^2*(1./l.^2 + 1i./(gp*l)));
for p = 1:2
    e = e + A(p)/li(p)*( exp(1i*ph(p))./(1/li(p) - 1./l - 1i/gi(p)) ...
                       + exp(-1i*ph(p))./(1/li(p) + 1./l + 1i/gi(p)) );
end
end
