function [ac_smooth aci] = smooth_ac(ac_prof,ac_alt,AK,prior,air,alt,varargin)
% [ac_smooth, aci] = smooth_ac(ac_prof,ac_alt,AK,prior,air,alt,remove_extras)
% if remove_extras = 1, values above and below original aircore are nan 

% extrapolate aircore
aci = extrapolate_ac(ac_prof,ac_alt,air,alt)./air*1e9; % ppb

prior = prior./air*1e9; % ppb

ac_smooth = prior(:) + AK*(aci(:)-prior(:)); % ppb

% remove extrapolated parts?
if (nargin==7 && varargin{1}==1)

    ind1 = find(ac_prof(2:end)>0,1,'first');
    ind2 = find(ac_prof>0,1,'last');

    ac_smooth(alt>ac_alt(ind1)) = nan;
    ac_smooth(alt<ac_alt(ind2)) = nan;

end