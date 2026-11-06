clear; 
set(0, 'DefaultFigureCreateFcn', @(h,~) set(h, 'Renderer', 'painters'))

%% Green's function check: Ewald vs direct spectral sum at |Y| >= 0.15 d
d = 1; k = 2*pi/0.7; 
beta = k*sin(20*pi/180);
X = rand(1,200)*2-1;  
Y = (rand(1,200)*0.6+0.15).*sign(rand(1,200)-0.5);
G = greenper1d(X, Y, k, beta, d, struct('regular',false));
Gs = zeros(size(X));
for n = -400:400
    bn = beta + 2*pi*n/d;  
    gn = sqrt(k^2 - bn^2); 
    if imag(gn) < 0 
        gn = -gn; 
    end
    Gs = Gs + (1i/(2*d))*exp(1i*bn*X + 1i*gn*abs(Y))/gn;
end
fprintf('Ewald vs spectral sum error: %.1e \n', max(abs(G-Gs))/max(abs(Gs)));