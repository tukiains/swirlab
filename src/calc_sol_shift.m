function shift = calc_sol_shift(wn,refe,sol,wn_shift,linewn)
% shift = calc_sol_shift(wn,refe,sol,wn_shift,linewn)
%
% finds solar wavelength shift by scanning

sol = conv_spectrum(wn,sol);

sol = interp1(wn,sol,wn+wn_shift,'linear','extrap');

r = refe./max(refe);

solind = findnearest(linewn,wn);

ind = find(wn>wn(solind)-0.06 & wn<wn(solind)+0.06);

dl = -0.03:0.0001:0.03;

for n=1:length(dl)
    sol2 = interp1(wn+dl(n),sol,wn(ind));
    diffu = sol2-r(ind);
    diffu = diffu.^2;
    diffz(n) = sqrt(mean(diffu));
end

ind = find(diffz==min(diffz));
shift = dl(ind);

