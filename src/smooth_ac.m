function [ac_smooth aci] = smooth_ac(ac_prof,ac_alt,AK,prior,air,alt,varargin)
% [ac_smooth, aci] = smooth_ac(ac_prof,ac_alt,AK,prior,air,alt,remove_extras)
% if remove_extras = 1, values above and below original aircore are nan 


pa = ac_alt;
p = ac_prof;

% extrapolate aircore
pa(1) = max(alt);
p(1) = 1e-12;
pa(2) = 35;
p(2) = 1e-12;

% also lower part 
pa(end) = 1e-12;
p(end) = p(find(p>0 & pa>0,1,'last'));

ind = find(p>0 & pa>0);
p = p(ind);
pa = pa(ind);

aci = interp1(pa,p,alt,'linear','extrap').*air;

% smooth 
ac_smooth = prior(:) + AK*(aci(:)-prior(:));

% remove extrapolated parts?
if (nargin==7 && varargin(1)==1)

    ind1 = find(ac_prof(2:end)>0,1,'first');
    ind2 = find(ac_prof>0,1,'last');

    ac_smooth(alt>ac_alt(ind1)) = nan;
    ac_smooth(alt<ac_alt(ind2)) = nan;

end