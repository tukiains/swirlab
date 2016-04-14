
clear all
close all

% set path for absorption coeffs:
voigt_path = '/home/tukiains/Dropbox/voigt_shapes/';

[pathstr,name] = fileparts(which('get_ftir_files.m'));

% all files in input folder
prefix = [pathstr,'/../input_data/ftir_spectra/'];
mfiles = dir([prefix,'so*']);
mfiles = struct2cell(mfiles);
mfiles = mfiles(1,:);

% what measurement file we study?
n = 3;
mfile = [prefix mfiles{n}]

% number of components
k = 3;

per_alt = 1:1:60;
per = 0.0001;

% reference
out = ch4_simu(voigt_path,mfile,k,per_alt(n),0);

for n = 1:length(per_alt)

    % perturbate:
    out_per(n) = ch4_simu(voigt_path,mfile,k,per_alt(n),per);
    
end
