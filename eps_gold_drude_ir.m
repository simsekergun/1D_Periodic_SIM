function e = eps_gold_drude_ir(lam_nm)
%EPS_GOLD_DRUDE_IR  Drude permittivity of gold in the near/mid infrared.
%   e = eps_gold_drude_ir(lam_nm),  exp(-i w t) convention (Im e > 0).
%   e = 1 - wp^2/(w (w + i g)),  w = 1e7/lam_nm (cm^-1),
%   wp = 7.28e4 cm^-1, g = 215 cm^-1 (Ordal et al., Appl. Opt. 24, 4493 (1985)).
%   Intended for lambda >~ 1 um.  Porto et al. (PRL 83, 2845, 1999) used
%   tabulated gold data instead, so small differences are expected.
w  = 1e7./double(lam_nm);
wp = 7.28e4;  g = 215;
e  = 1 - wp^2./(w.*(w + 1i*g));
end
