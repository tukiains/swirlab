function shift = calc_sol_shift(mfile,wn_shift,labpath,solar_line_file,sol_shift_wn)

% measurement
[t,wn] = read_ftir_spectrum(mfile,[sol_shift_wn-5 sol_shift_wn+5]);
f = ggg_ils(2,length(wn),median(wn),median(diff(wn)),45,0);
t = conv(t,f,'same')./max(t);
t = interp1(wn+wn_shift,t,wn,'linear','extrap');

% sol
sol = calc_solar_spectrum(wn,[labpath,'/../',solar_line_file]);
f = ggg_ils(2,length(wn),median(wn),median(diff(wn)),45,2.3923e-3);
sol = conv(sol,f,'same');

% only the strong peak
ind = find(wn>sol_shift_wn-0.15 & wn<sol_shift_wn+0.15);
t = t(ind);
t = t/max(t);
wn = wn(ind);

sol = sol(ind);
sol = sol/max(sol);

dl = -0.03:0.0001:0.03;

for n = 1:length(dl)
    sol2 = interp1(wn+dl(n),sol,wn,'linear','extrap');
    ss(n) = sum((sol2-t).^2);
end

ind = find(ss==min(ss));
shift = dl(ind);

