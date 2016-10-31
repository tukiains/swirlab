function [voi]=voigt_shape(wl,I,grid,bP,bD)
%
% voi = voigt_shape(wl,I,grid,bP,bD)
%
% INPUT:
% wl   = wavenumber of the line (1/cm)
% I    = intensity of the line
% grid = target grid (1/cm)
% dP   = Pressure broadening half-width. Lorentzian shape. Pressure.
% bD   = Doppler broadening half-width. Gaussian shape. Temperature.
%
% Voigt profile approximation taken from:
% Gharavi, Mohammadreza; Buckley, Steven.
% Single Diode Laser Sensor for Wide-Range H2O Temperature Measurements.
% Applied spectroscopy, V58,468(2004).

gV  = 0.5346*bP + sqrt(0.2166*bP^2 + bD^2); % Voigt profile half width
x   = bP/gV;
y   = abs(grid-wl)/gV; 
S1  = I/(2*gV*(1.065 + 0.447*x + 0.058*x^2));
voi = S1*((1-x)*exp(-0.693.*y.^2) + (x./(1+y.^2)) + ...
          0.016*(1-x)*x*(exp(-0.0841.*y.^2.25)-1./(1 + 0.021.*y.^2.25)));
    
