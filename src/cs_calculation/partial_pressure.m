function pp = partial_pressure(p,ngas,nair,varargin)
% pp = partial_pressure(p,ngas,nair,unit)
%
% INPUT:
% p = pressure [hPa]
% ngas = number density of gas
% nair = number density of air
% unit can be 'atm' or 'pascal' or empty when it's hPa

p = p(:);
ngas = ngas(:);
nair = nair(:);

pp = ngas./nair.*p;

if (nargin==3)
    return;
else
    if (varargin{1}=='atm')
        pp = pp/1013.25;
    elseif (varargin{1}=='pascal')
        pp = pp*1000;
    else
        error('Unknown unit')
    end
end