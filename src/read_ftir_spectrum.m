function [refe,refe_wn] = read_ftir_spectrum(fname,wnrange);
% [refe,refe_wl] = read_ftir_spectrum(fname,wnrange)

a = textread(fname);
refe_wn = a(:,1); 
refe = a(:,2);

ind = find(refe_wn>min(wnrange) & refe_wn<max(wnrange));

refe = refe(ind);
refe_wn = refe_wn(ind);

