function [sza names year day hour lat lon alt tins pins hins tout pout hout sia] = read_zenith_angles(fname,zenlim)
%
% read GGG2014 .grl file
    
% number of header lines
fid = fopen(fname);

nheader = str2num(fgetl(fid));
nheader = nheader(1);

% read file
s = textscan(fid, ...
             ['%s ' ... % Spectrum_File_Name
              '%d32 ' ... % Year
              '%d32 ' ... % Day
              '%f32 ' ... % Hour
              '%f32 ' ... % oblat
              '%f32 ' ... % oblon
              '%f32 ' ... % obalt
              '%f32 ' ... % ASZA   
              '%f32 ' ... % POFF
              '%f32 ' ... % AZIM
              '%f32 ' ... % OSDS
              '%f32 ' ... % OPD
              '%f32 ' ... % FOVI
              '%f32 ' ... % FOVO
              '%f32 ' ... % AMAL
              '%d32 ' ... % IFIRST
              '%d32 ' ... % ILAST
              '%f32 ' ... % DELTA_NU
              '%d32 ' ... % POINTER
              '%d32 ' ... % BPW 
              '%f32 ' ... % ZOFF
              '%d32 ' ... % SNR
              '%s '   ... % APF
              '%f32 ' ... % tins
              '%f32 ' ... % pins
              '%f32 ' ... % hins
              '%f32 ' ... % tout
              '%f32 ' ... % pout
              '%f32 ' ... % hout
              '%f32 ' ... % sia
              '%f32 ' ... % fvsi
              '%f32 ' ... % wspd
              '%f32 ' ... % wdir
              '%f32 ' ... % lasf
              '%f32 ' ... % wavtkr
              '%f32'], ... % aipl
             'headerlines',nheader-1); 

fclose(fid);

% screen by sza
sza = s{8};
ind = find(sza<zenlim);
sza = sza(ind);

names = cell2mat(s{1}(ind));
year = s{2}(ind);
day = s{3}(ind);
hour = s{4}(ind);
lat = s{5}(ind);
lon = s{6}(ind);
alt = s{7}(ind);
azi = s{10}(ind);
tins = s{24}(ind);
pins = s{25}(ind);
hins = s{26}(ind);
tout = s{27}(ind);
pout = s{28}(ind);
hout = s{29}(ind);
sia = s{30}(ind);

