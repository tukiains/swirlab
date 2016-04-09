
clear all
close all

% set path for absorption coeffs:
voigt_path = '/home/tukiains/Dropbox/voigt_shapes/';

[pathstr,name] = fileparts(which('get_ftir_files.m'));

lm_only = false;
lis = false;
k = 3;

% all files we have
prefix = [pathstr,'/../input_data/ftir_spectra/'];
mfiles = dir([prefix,'so*','*.0*']);

for n = 1:length(mfiles)
    mfile(n,:) = [prefix,mfiles(n).name];
    out(n) = ftir_dimred_mcmc(voigt_path,mfile(n,:),lm_only,lis,k);
end

%save('/home/tukiains/Dropbox/Public/ftir_results.mat','out','mfile');