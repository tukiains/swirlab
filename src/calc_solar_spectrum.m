function sol = calc_solar_spectrum(wn,fname)
% sol = calc_solar_spectrum(wn,solar_line_file)

a     = importdata(fname);
sol   = interp1(a(:,1),a(:,2),wn);
M     = planck(5778,wn2wl(wn)); 
sol   = sol(:).*M(:);



