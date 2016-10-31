function bP = calc_pressure_broadening(T,P,Ps,n,s_air,s_self);
% [bP] = calc_pressure_broadening(T,P,Ps,n,s_air,s_self);
%
% INPUT:
% T = temperature [K]
% P = pressure [atm]
% Ps = partial pressure [atm]
% n = temperature dependency coeff.
% s_air = air broaneded halfwidth
% s_self = self broaneded halfwidth 
%
% OUTPUT:
% pP = Pressure broadening
%
% NOTES:
% reference T = 296 K and reference P = 1 atm (HITRAN)

bP = (296/T).^n.*(s_air.*(P-Ps)+s_self.*Ps); % Pressure broadening
