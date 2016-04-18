
clear all
close all

% set path for absorption coeffs:
voigt_path = '/home/tukiains/Dropbox/voigt_shapes/';

[pathstr,name] = fileparts(which('get_ftir_files.m'));

% all files in input folder
%prefix = [pathstr,'/../input_data/ftir_spectra/'];
prefix = '/home/tukiains/Documents/ggg-2014/ggg-stable/i2s/opus-i2s/spectra/20140319/';
mfiles = dir([prefix,'so*']);
mfiles = struct2cell(mfiles);
mfiles = mfiles(1,:);

% what measurement file we study?
n = 50;
mfile = [prefix mfiles{n}]

% number of components
k = 3;

per_alt = 1:10:70;
per = 1e10; % 1 % perturbation

% starting point
out = ch4_simu(voigt_path,mfile,k,0,0);

refe = out.dr_lm_atmos.ch4;

refe_col = trapz(out.geo.center_alts,refe);

refe_col_true = trapz(out.simuatmos.alt,out.simuatmos.ch4);

for n = 1:length(per_alt)

    % perturbate:
    out_per = ch4_simu(voigt_path,mfile,k,per_alt(n),per);
    
    col = trapz(out_per.geo.center_alts,out_per.dr_lm_atmos.ch4);

    A(n,:) = (out_per.dr_lm_atmos.ch4-refe)./out_per.simuatmos.dx;
    
    %Ac(n) = (col-refe_col)./out_per.simuatmos.dx;
    
    Ac(n) = (col-refe_col)/(trapz(out_per.simuatmos.alt,out_per.simuatmos.ch4)-refe_col_true);

end


figure(1)
plot(Ac,per_alt,'ro-')
hold on
plot(out.A_column,out.geo.center_alts,'bo-');
legend('Numerical','Analytic')

% true column

%refe_col_dr = out.dr_
%col_refe = trapz(out
