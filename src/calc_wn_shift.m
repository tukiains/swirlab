function shift = calc_wn_shift(geo,wn,gasvec,cros,refe,sol,L,ncut)
% calc_wn_shift(geo,wn,gasvec,cros,refe,sol,L,ncut)
%
% finds wavelength shift between data and model by simple scanning

% run forward model
t = calc_direct_radiance(geo.layer_dens,geo.los_lens,gasvec,cros,sol,1,1,1,0,L); 

% convolute
tc = conv_spectrum(wn,t);
tc = tc(ncut:end-ncut);
refe = refe(ncut:end-ncut);
wn = wn(ncut:end-ncut);

% approximate fit
tc2 = tc*(tc\refe);

% investigate only 50 smallest points 
[~,ind] = sort(tc2);
ind = ind(1:50);

tc3 = tc2(ind);
refe2 = refe(ind);
wn2 = wn(ind);

dl = -0.01:0.00001:0.015;

for n=1:length(dl)
    tc4 = interp1(wn2,tc3,wn2+dl(n));
    ind = find(isnan(tc4)==0);
    if (ind>0)
        diffu = tc4(ind)-refe2(ind);
        diffu = diffu.^2;
        diff(n) = sqrt(mean(diffu));
    else
        diff(n) = 1;
    end
end

ind = find(diff==min(diff));
shift = dl(ind);

