
clear all
close all

% set path for absorption coeffs:
voigt_path = '/home/tukiains/data/voigt_shapes/';

[pathstr,name] = fileparts(which('get_ftir_files.m'));

% two test files
mfile = [pathstr, '/../input_data/ftir_spectra/wg20140319saebaa.067'];
%mfile = [pathstr, '/../input_data/ftir_spectra/wg20130903saebaa.127'];

% use likelihood-informed retrieval?
lis = true;

% retrieve ch4
out = ftir_dimred_mcmc(voigt_path,mfile,lis);

figure(1)
clf
hold on
fa = 0.4;

% prior
h1 = show_prior(out.geo.center_alts,out.geo.layer_dens.ch4,out.geo.air,out.dr_lm_P{1})
set(h1,'facealpha',fa)

% dimension reduction with MCMC (2-sigma posterior):
h2 = plot_curtain(out.geo.center_alts,plims(out.mcmc_profs,[0.025 0.5 0.975]),[.5 .7 .3]);
set(h2,'facealpha',fa)

% scaling of prior profile:
%plot(out.geo.layer_dens.ch4./out.geo.air*out.scaling_factors(1),out.geo.center_alts,'b-','linewidth',2)

% dimension reduction with LM:
h4 = show_lm(out.geo.center_alts,out.geo.layer_dens.ch4,out.geo.air,out.dr_lm_P{1},out.dr_lm_theta,out.dr_lm_cmat,out.dr_k);
set(h4,'facealpha',fa)

% aircore
ac_file = get_aircore_file(mfile(end-17:end-10),[pathstr, '/../input_data/aircore/']);
if (exist(ac_file)==2)
    [co2,co2e,ch4,ch4e,co,coe,pres,alt,temp,air] = read_aircore_sounding(ac_file);
    plot(ch4./1e9,alt,'r-','linewidth',2)
end

set(gca,'ylim',[0 45])

l = legend('prior','DR (MCMC)', 'DR (LM)','AirCore')
legend boxoff



