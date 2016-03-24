function [sza names year lat lon alt tins pins hins tout pout hout sia] = read_zenith_angles(fname,zenlim)
% fname = *.grl file

fid = fopen(fname);

s = textscan(fid, ...
             ['%s %d32 %d32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %d32 %d32 %f32 %f32 ' ...
              '%d32 %f32 %d32 %s %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32 %f32'], ...
             'headerlines',3); 

fid = fclose(fid);

names = cell2mat(s{1});
year = s{2};
day = s{3};
hour = s{4};
lat = s{5};
lon = s{6};
alt = s{7};
sza = s{8};
azi = s{10};
tins = s{24};
pins = s{25};
hins = s{26};
tout = s{27};
pout = s{28};
hout = s{29};
sia = s{30};

ind = find(sza<zenlim);
names = names(ind,:);
year = year(ind);
day = day(ind);
hour = hour(ind);
lat = lat(ind);
lon = lon(ind);
alt = alt(ind);
sza = sza(ind);
azi = azi(ind);
tins = tins(ind);
pins = pins(ind);
hins = hins(ind);
tout = tout(ind);
pout = pout(ind);
hout = hout(ind);
sia = sia(ind);