function x = pressure_shift(xin,sigma,p)
% x = pressure_shift(xin,sigma,p)
%
% correction for the line position due to pressure
%
% INPUT:
%
% xin = wavenumber (1/cm)
% sigma = Air broadened pressure shift
% p = pressure [atm]
%
% OUTPUT:
%
% x = corrected line position

x = xin + sigma.*p;