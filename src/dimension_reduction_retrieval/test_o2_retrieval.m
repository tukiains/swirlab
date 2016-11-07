
clear all
close all

% set path for absorption coeffs:
voigt_path = '/home/tukiains/data/voigt_shapes/';

[pathstr,name] = fileparts(which('get_ftir_files.m'));

% all files in input folder
prefix = [pathstr,'/../input_data/ftir_spectra/'];

mfiles = dir([prefix,'so*']);
mfiles = struct2cell(mfiles);
mfiles = mfiles(1,:);

n = 1;
mfile = [prefix mfiles{n}]

zenlim = 82;

% retrieve o2
out = o2_retrieval(voigt_path,mfile,zenlim);

figure(1)
clf
hold on

% prior
plot(out.geo.layer_dens.o2./out.geo.air,out.geo.center_alts,'b--','linewidth',2)

% scaled prior
plot(out.geo.layer_dens.o2./out.geo.air*out.scaling_factors(1),out.geo.center_alts,'b-','linewidth',2)

figure(2)
clf
plot(out.wn,out.scaling_residual,'r-')

