function sza_out = get_sza_angle(mfile,grlfile)
% sza = get_sza_angle(mfile,grlfile)
% 
% read solar zenith angle of the ftir measurement 
% needs the grl file, too

[sza names year lat lon alt tins pins hins tout pout hout sia] = ...
    read_zenith_angles(grlfile,90);

runs = str2num(names(:,end-3:end));

runno = str2num(mfile(end-2:end));

ind = find(runs==runno);

sza_out = double(sza(ind(1)));