function shift = calc_wn_shift(geo,wn,gasvec,cros,refe,sol,mindep,L)
% calc_wn_shift(geo,wn,gasvec,cros,refe,sol,mindep,L)
%
% finds wavelength shift between data and model by simple scanning

% run forward model
[t,~,~] = calc_direct_radiance(geo,wn,gasvec,cros,sol,1,1,1,0,L); 

% convolute
tc = conv_spectrum(wn,t);

% approximate fit
tc2 = tc*(tc\refe);

ind = find(tc2<mindep); % investigate only peaks 
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

