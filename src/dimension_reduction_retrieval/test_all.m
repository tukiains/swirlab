
clear all
close all

% set path for absorption coeffs:
voigt_path = '/home/tukiains/Dropbox/voigt_shapes/';

[pathstr,name] = fileparts(which('get_ftir_files.m'));

usesimu = false;
fixo = true; % fix offset in prior scaling?
lm_only = true;
lis = false;
k = 3;
mdate = '20140715';

% all files we have
%prefix = [pathstr,'/../input_data/ftir_spectra/'];
prefix = ['/home/tukiains/Documents/ggg-2014/ggg-stable/i2s/opus-i2s/spectra/',mdate,'/'];

mfiles = dir([prefix,'so*','*.0*']);

g = 1;
for n = 1:length(mfiles)
    fname = [prefix,mfiles(n).name];
    temp = ftir_dimred_mcmc(voigt_path,fname,lm_only,lis,k,fixo,usesimu);
    if (temp.sza>0)
        mfile(g,:) = fname;
        out(g) = temp;
        g = g+1;
    end
end

save(['/home/tukiains/Dropbox/Public/ftir_results_',mdate,'_fixo.mat'],'out','mfile');