
clear all
close all

% set path for absorption coeffs:
voigt_path = '/home/tukiains/Dropbox/voigt_shapes/';

% all files in input folder
mdate = '20130903';
prefix = ['/home/tukiains/Documents/ggg-2014/ggg-stable/i2s/opus-i2s/spectra/',mdate,'/'];

files = list_files(prefix,'so*.0*');

% zenith angle limit
zenlim = 90;
% how many principal components
k = 3;
% fixed offset?
fixo.scale = true;
fixo.dr = true;

out = retrieve_all(voigt_path,files,k,fixo,zenlim);


