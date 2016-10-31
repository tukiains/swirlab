function [bD] = calc_doppler_broadening(T,freq,M);
% [bD] = calc_doppler_broadening(T,freq,M);
%
% INPUT
% T = temperature [K]
% freq = wavenumber of the line [1/cm]
% M = molecular weight
%
% OUTPUT
% bD = Doppler broadening

bD = 3.581e-7.*freq*sqrt(T/M); 

