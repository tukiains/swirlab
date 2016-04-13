function ac_smooth = smooth_ac(ac_prof,ac_alt,AK,prior,air,alt)

% extrapolate aircore
ac_prof(1) = 1e-12;
ac_alt(1) = max(alt);
ind = find(ac_prof>0 & ac_alt>0);
aci = interp1(ac_alt(ind),ac_prof(ind),alt,'linear','extrap').*air;

% smooth 
ac_smooth = prior(:) + AK*(aci(:)-prior(:));

ind1 = find(ac_prof(2:end)>0,1,'first');
ind2 = find(ac_prof>0,1,'last');

ac_smooth(alt>ac_alt(ind1)) = nan;
ac_smooth(alt<ac_alt(ind2)) = nan;
