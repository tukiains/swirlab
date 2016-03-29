function [refe,refe_wn] = read_ftir_spectrum(fname,wnrange);
% [refe,refe_wl] = read_ftir_spectrum(fname,wnrange)

fid = fopen(fname);
data = textscan(fid, '%f64%f64');

refe_wn = data{1}; 
refe = data{2};

ind = find(refe_wn>min(wnrange) & refe_wn<max(wnrange));

refe = refe(ind);
refe_wn = refe_wn(ind);

fclose(fid);