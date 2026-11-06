function w = faddeeva_w(z, N)
%FADDEEVA_W  Faddeeva function w(z) = exp(-z^2) erfc(-i z) for complex z.
%   Weideman's rational (FFT) approximation, J.A.C. Weideman, SIAM J. Numer.
%   Anal. 31, 1497 (1994), with reflection w(z) = 2 exp(-z^2) - w(-z) for
%   Im z < 0.  N = 40 (default) gives ~3e-14 relative accuracy for |z| < 7.
%   Then erfc(z) = exp(-z^2) w(i z).
if nargin < 2, N = 40; end
persistent aN LN NN
if isempty(NN) || NN ~= N
    M = 2*N;  
    M2 = 2*M;  
    k = (-M+1:M-1).';
    L = sqrt(N/sqrt(2));
    t = L*tan(k*pi/M/2);
    f = exp(-t.^2).*(L^2 + t.^2);  
    f = [0; f];
    a = real(fft(fftshift(f)))/M2;
    aN = flipud(a(2:N+1));  
    LN = L;  
    NN = N;
end
neg = imag(z) < 0;
zz = z;  
zz(neg) = -z(neg);
Z = (LN + 1i*zz)./(LN - 1i*zz);
p = polyval(aN, Z);
w = 2*p./(LN - 1i*zz).^2 + (1/sqrt(pi))./(LN - 1i*zz);
w(neg) = 2*exp(-z(neg).^2) - w(neg);
end
