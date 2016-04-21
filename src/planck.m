function [I_out] = planck(T,wl)
%
% I = planck(T,wl)
%
% wl = [nm]
% I  = [W/sr/m^2/nm]

wl = wl*1e-9; % nm -> m

k = 1.38064853e-23;
h = 6.626070041e-34;
c = 299792458;

for n=1:length(wl)
    a = 2*h*c^2/wl(n)^5;
    b = h*c / (wl(n)*k*T);
    I_out(n) = a * (1/(exp(b)-1)); % W/sr/m^2/m
end

I_out = I_out*1e-9; % to units we want


