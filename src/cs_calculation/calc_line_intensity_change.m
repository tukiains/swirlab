function [line_out] = calc_line_intensity_change(freq,line_in,El,Q_ref,Qt,T)
%  [line_out] = calc_line_intensity_change(freq,line_in,El,Q_ref,Qt,T)

% freq = line frequency [1/cm]
% line_in = line intensity in T_ref (296 K)
% El = lower state energy of the transition
% Q_ref = Q value in T_ref
% Qt = Q value in T
% T = temperature

T_ref = 296; % [K]
c2 = 1.4388; % [cm K]

F1 = Q_ref/Qt;

F2 = (exp(-c2.*El/T))./(exp(-c2.*El/T_ref));

F3 = (1-exp(-c2.*freq/T))./(1-exp(-c2.*freq/T_ref));

line_out = line_in.*F1.*F2.*F3;