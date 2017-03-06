
clear all
close all

% set path for absorption coeffs:
voigt_path = '/home/tukiains/Dropbox/voigt_shapes/';

% all files in input folder
mdate = '20140319';
%prefix = ['/home/tukiains/Documents/ggg-2014/ggg-stable/i2s/opus-i2s/spectra/',mdate,'/'];
prefix = ['/media/tukiains/3b48ca75-37ff-42ba-ab71-b78f39cd9a79/tukiains/Documents/ggg-2014/ggg-stable/i2s/opus-i2s/spectra/',mdate,'/'];

files = list_files(prefix,'so*.0*');

% zenith angle limit
zenlim = 90;
% how many principal components
k = 3;
% fixed offset?
fixo.scale = true;
fixo.dr = true;

day = retrieve_all(voigt_path,files,k,fixo,zenlim);

figure(1)
clf
hold on
for n=1:length(day.profiles)
    plot(day.profiles(n).dr_lm_atmos.ch4./day.air,day.center_alts,'g-')
end

% aircore
[pathstr,name] = fileparts(which('get_ftir_files.m'));
ac_file = get_aircore_file(mdate,[pathstr,'/../input_data/aircore/']);

[co2,co2e,ch4,ch4e,co,coe,pres,alt,temp,air] = read_aircore_sounding(ac_file);
ch4 = ch4/1e9;
plot(ch4,alt,'r-','linewidth',2)


