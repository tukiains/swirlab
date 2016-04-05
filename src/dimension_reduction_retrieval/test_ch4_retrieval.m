
clear all
close all

% path for absorption coeffs:
voigt_path = '/foo/bar/voigt_shapes/';

% measurement file:
%mfile = '/foo/bar/swirlab/input_data/ftir_spectra/wg20140319saebaa.067';
mfile = '/foo/bar/swirlab/input_data/ftir_spectra/wg20130903saebaa.127'

% use likelihood-informed retrieval?
lis = false;

% retrieve ch4
out = ftir_dimred_mcmc(voigt_path,mfile,lis);

figure(1)
clf
hold on
% dimension reduction with MCMC (2-sigma posterior):
plot_curtain(out.geo.center_alts,plims(out.mcmc_profs,[0.025 0.5 0.975]),[.5 .7 .3],[.5 .7 .3]);
% scaling of prior profile:
plot(out.geo.layer_dens.ch4./out.geo.air*out.scaling_factors(1),out.geo.center_alts,'b-','linewidth',2)
% dimension reduction with LM:
plot(out.dr_lm_atmos.ch4./out.geo.air,out.geo.center_alts,'g-','linewidth',2)
set(gca,'ylim',[0 40])

