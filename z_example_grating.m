% Gold elliptical nanowire grating (TE_z), R/T/A spectrum.
clear;


geom = struct('type','ellipse','a',80,'b',40);    % nm
d    = 400;                                         % period (nm)
lams = 450:5:900;
R = zeros(size(lams)); 
T = R; 
A = R; 
res = R;
for k = 1:numel(lams)
    o = sim1dperiodic(geom, lams(k), eps_gold(lams(k)), 96, d, struct('pol','TE','thetai',10*pi/180));
    R(k) = o.R; 
    T(k) = o.T; 
    A(k) = o.A; 
    res(k) = o.energy_residual;
end
fprintf('max |1-R-T-A| = %.1e\n', max(abs(res)));

figure(21); clf;
plot(lams, R, lams, T, lams, A); 
grid on;
axis tight;
legend('R','T','A', 'Location', 'East'); 
xlabel('\lambda (nm)');
title('Gold elliptical wires, d = 400 nm, \theta_i = 10^\circ, TE_z');
