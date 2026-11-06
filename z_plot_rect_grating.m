% Three periods of a grating of nearly rectangular cylinders.
%   Uses superellipse_geom.m 
%   Left: three unit cells for p = 16, with the period d and the unit-cell
%   boundaries |x| = d/2.  Right: one cylinder for p = 8, 16, 32.

clear; clc;
a = 100;        % half-width
b = 50;          % half-height (nm)
d = 400;        % period (nm), must satisfy 2a < d
p = 16;          % superellipse exponent (even)
t = linspace(0, 2*pi, 721);

fig = figure('Color','w','Position',[100 100 900 330]);

% ---- (a) three periods ------------------------------------------------
ax1 = subplot(1,3,[1 2]); 
hold(ax1,'on');
geom = superellipse_geom(a, b, p);
[r, ~, ~] = geom.fun(t);
gold = [0.83 0.69 0.22];
for m = -1:1
    fill(ax1, r(1,:) + m*d, r(2,:), gold, 'EdgeColor','k', 'LineWidth',1);
    xline(ax1, m*d - d/2, ':', 'Color',[0.5 0.5 0.5]);
end
xline(ax1, 1.5*d, ':', 'Color',[0.5 0.5 0.5]);
% period annotation
plot(ax1, [0 d], -1.6*b*[1 1], 'k-', 'LineWidth',0.8);
plot(ax1, [0 0; d d].', [-1.75*b -1.45*b; -1.75*b -1.45*b].', 'k-');
text(ax1, d/2, -1.95*b, sprintf('d = %d nm', d), 'HorizontalAlignment','center');
% incident wave from above (theta_i = 0)
quiver(ax1, -d, 2.6*b, 0, -0.9*b, 0, 'k', 'LineWidth',1.2, 'MaxHeadSize',0.8);
text(ax1, -d + 12, 2.4*b, 'incident', 'HorizontalAlignment','left');
axis(ax1, 'equal'); 
xlim(ax1, [-1.6*d 1.6*d]); 
ylim(ax1, [-2.9*b 3*b]);
xlabel(ax1, 'x (nm)'); 
ylabel(ax1, 'y (nm)');
title(ax1, sprintf('(a) three periods, %d \\times %d nm, p = %d', 2*a, 2*b, p));
box(ax1,'on');

% ---- (b) effect of p ---------------------------------------------------
ax2 = subplot(1,3,3); 
hold(ax2,'on');
ps = [8 16 32];  
cols = lines(numel(ps));
for k = 1:numel(ps)
    g = superellipse_geom(a, b, ps(k));
    [rk, ~, ~] = g.fun(t);
    plot(ax2, rk(1,:), rk(2,:), 'Color', cols(k,:), 'LineWidth', 1.2);
end
plot(ax2, [-a a a -a -a], [-b -b b b -b], 'k--');      % ideal rectangle
legend(ax2, [compose('p = %d', ps), {'rectangle'}], 'Location','southoutside', ...
       'NumColumns', 3, 'Box','off');
axis(ax2, 'equal'); 
xlim(ax2, [-1.15*a 1.15*a]);
xlabel(ax2, 'x (nm)'); 
title(ax2, '(b) cross-section vs. p'); 
box(ax2,'on');
% exportgraphics(fig, 'figure_rect_grating.png', 'Resolution', 300);
