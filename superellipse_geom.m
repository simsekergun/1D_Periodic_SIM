function geom = superellipse_geom(a, b, p)
% Analytic, nearly rectangular boundary for sim1dperiodic / sim2d.
%   geom = superellipse_geom(a, b, p)
%
%   Polar superellipse   rho(t) = ( cos(t)^p/a^p + sin(t)^p/b^p )^(-1/p),
%   r(t) = rho(t) [cos t; sin t],  t in [0, 2*pi), counter-clockwise.
%   a, b : half-width (x) and half-height (y) of the bounding rectangle
%   p    : EVEN integer (4, 8, 16, 32, ...).  For even p the curve is
%          real-analytic, so SIM keeps exponential convergence; p -> Inf
%          gives the 2a x 2b rectangle.  The corner radius of curvature
%          scales as 1/p (about 187/p nm for a = 100, b = 50 nm).
%
%   Returns geom.type = 'custom' and geom.fun = @(t) [r, dr, ddr] (2xN each),
%   with analytic first and second derivatives (verified against finite
%   differences to the truncation level).
%
%   Periodic use: the cross-section must lie inside the strip |x| < d/2,
%   i.e. 2a < d.
%
% Ergun Simsek, UMBC

if mod(p,2) ~= 0 || p < 2
    error('superellipse_geom:p','p must be an even integer >= 2');
end
geom = struct('type','custom','a',a,'b',b,'p',p);
geom.fun = @(t) curve_superellipse(t, a, b, p);
end

function [r, dr, ddr] = curve_superellipse(t, a, b, p)
t  = t(:).';
c  = cos(t);  
s = sin(t);
q  = c.^p/a^p + s.^p/b^p;
q1 = p*( -c.^(p-1).*s/a^p + s.^(p-1).*c/b^p );
q2 = p*( (p-1)*c.^(p-2).*s.^2/a^p - c.^p/a^p ...
       + (p-1)*s.^(p-2).*c.^2/b^p - s.^p/b^p );
rho = q.^(-1/p);
r1  = -(1/p)*q.^(-1/p-1).*q1;
r2  =  (1/p)*(1/p+1)*q.^(-1/p-2).*q1.^2 - (1/p)*q.^(-1/p-1).*q2;
r   = [rho.*c;               rho.*s];
dr  = [r1.*c - rho.*s;       r1.*s + rho.*c];
ddr = [(r2 - rho).*c - 2*r1.*s;  (r2 - rho).*s + 2*r1.*c];
end
