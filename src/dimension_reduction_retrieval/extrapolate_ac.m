function aci = extrapolate_ac(ac_prof,ac_alt,air,alt)
% aci = extrapolate_ac(ac_prof,ac_alt,alt)
%
% ac_prof: aircore profile [mole fraction]
% ac_alt: aircore altitude [km]
% air: density
% alt: target altitude

pa = ac_alt;
p = ac_prof;

% extrapolate aircore
pa(1) = max(alt);
p(1) = 1e-12;
pa(2) = 60;
p(2) = 1e-12;

% also lower part
pa(end) = 1e-12;
p(end) = p(find(p>0 & pa>0,1,'last'));

ind = find(p>0 & pa>0);
p = p(ind);
pa = pa(ind);

aci = interp1(pa,p,alt,'linear','extrap').*air;
