%EXAMPLE_RECTANGLE  Nearly rectangular gold grating with sim1dperiodic.
%   Convergence in N for several superellipse exponents p.
%   Choose N with the curvature term N_g = c_g*|dOmega|*max|kappa| (c_g = 2.5)
%   as a starting point, then confirm with a doubling test.
%   Note: the convergence rate slows as p grows, because the curve's
%   analyticity strip narrows (corner radius ~ 1/p).
%
%   Expected (NumPy transcription, gold 633 nm, 200x100 nm, d = 400 nm, TE,
%   normal incidence; |dR|, |dT| against N = 768):
%     p =  4: R = 0.498230827109, T = 0.470294172780;  <1e-12 at N = 128
%     p =  8: R = 0.529751985130, T = 0.436070211502;  <1e-12 at N = 256
%     p = 16: R = 0.543251605218, T = 0.420834715190;  3e-11 at N = 384
%     p = 32: R = 0.548921309160, T = 0.414184767253;  3e-8 at N = 512
%             (p = 32 needs N of about 1024 for 1e-10; use N = 1536 as reference)
%   The curvature estimate N_g (57, 120, 248, 506) underestimates the N
%   needed for 1e-10 by about a factor of 2 for these shapes.

clear; clc;
lam  = 633;  epsr = eps_gold(lam);       % make sure eps_gold.m is on the path
a = 100; b = 50;                         % 200 nm x 100 nm cross-section
d = 400;                                 % period (must satisfy 2a < d)
opts = struct('pol','TE','thetai',0);    % as for your other sim1dperiodic runs
ps = [4 8 16 32];
Ns = [64 96 128 192 256 384 512 768];

for p = ps
    geom = superellipse_geom(a, b, p);
    % geometric resolution estimate
    tt = linspace(0, 2*pi, 20001); tt(end) = [];
    [~, dr, ddr] = geom.fun(tt);
    J   = hypot(dr(1,:), dr(2,:));
    kap = (dr(1,:).*ddr(2,:) - dr(2,:).*ddr(1,:))./J.^3;
    L   = mean(J)*2*pi;
    fprintf('p = %2d: perimeter %.1f nm, min radius of curvature %.2f nm, N_g = %d\n', ...
            p, L, 1/max(kap), ceil(2.5*L*max(kap)));
    out = cell(numel(Ns),1);
    for k = 1:numel(Ns)
        out{k} = sim1dperiodic(geom, lam, epsr, Ns(k), d, opts);
    end
    % convergence of reflectance/transmittance against the largest N
    for k = 1:numel(Ns)-1
        fprintf('   N = %4d: |dR| = %.2e  |dT| = %.2e\n', Ns(k), ...
            abs(out{k}.R - out{end}.R), abs(out{k}.T - out{end}.T));
    end
end
