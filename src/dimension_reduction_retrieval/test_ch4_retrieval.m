
clear all
close all

% set path for absorption coeffs:
voigt_path = '/home/tukiains/data/voigt_shapes_2016/';

[pathstr,name] = fileparts(which('get_ftir_files.m'));

% all files in input folder
%prefix = [pathstr,'/../input_data/ftir_spectra/'];

prefix = '/home/tukiains/Documents/ggg-2014/ggg-stable/i2s/opus-i2s/spectra/20130903/';

mfiles = dir([prefix,'so*']);
mfiles = struct2cell(mfiles);
mfiles = mfiles(1,:);

n = 1;
mfile = [prefix mfiles{n}]

zenlim = 82;
usesimu = false;
lm_only = true;
lis = false;
k = 3;
% fixed offset?
fixo.scale = false; 
fixo.dr = false;
fixo.mcmc = false;

% retrieve ch4
out = ftir_dimred_mcmc(voigt_path,mfile,lm_only,lis,k,fixo,usesimu,zenlim);

figure(1)
clf
hold on
fa = 0.6;

% prior
h1 = show_prior(out.geo.center_alts,out.geo.layer_dens.ch4,out.geo.air,out.dr_lm_P{1})
set(h1,'facealpha',fa-0.2)

if (lm_only)
    % dimension reduction with LM:
    h2 = show_lm(out.geo.center_alts,out.geo.layer_dens.ch4,out.geo.air,out.dr_lm_P{1},out.dr_lm_theta,out.dr_lm_cmat,out.dr_k);
    set(h2,'facealpha',fa+0.2)
    str = 'DR (LM)';
else
    % dimension reduction with MCMC (2-sigma posterior):
    h2 = plot_curtain(out.geo.center_alts,plims(out.mcmc_profs,[0.025 0.5 0.975]),[.5 .7 .3]);
    set(h2,'facealpha',fa+0.2)
    if (lis)
        str = 'LIS (MCMC)';
    else
        str = 'DR (MCMC)';
    end
end

% prior
plot(out.geo.layer_dens.ch4./out.geo.air,out.geo.center_alts,'b--','linewidth',2)

% scaled prior
%plot(out.geo.layer_dens.ch4./out.geo.air*out.scaling_factors(1),out.geo.center_alts,'b-','linewidth',2)

% aircore
ac_file = get_aircore_file(get_date(mfile),[pathstr, '/../input_data/aircore/']);

[co2,co2e,ch4,ch4e,co,coe,pres,alt,temp,air] = read_aircore_sounding(ac_file);
ch4 = ch4/1e9;
plot(ch4,alt,'r-','linewidth',2)

% smoothed aircore
[ch4_smooth aci] = smooth_ac(ch4,alt,out.A_layer,out.geo.layer_dens.ch4,out.geo.air,out.geo.center_alts);
plot(ch4_smooth./out.geo.air',out.geo.center_alts,'k-','linewidth',2)

set(gca,'ylim',[0 30])
set(gca,'xlim',[.2e-6 2.1e-6])

%l = legend('prior',str,'prior mean','scaled prior','AirCore')
%legend boxoff

figure(2)
hold on
plot(aci./out.geo.air,out.geo.center_alts,'k-','linewidth',2)
plot(ch4,alt,'ro')


figure(3)
hold on
plot(out.geo.layer_dens.ch4./out.geo.air*out.scaling_factors(1),out.geo.center_alts,'b-','linewidth',2)
[ch4_smooth aci] = smooth_ac(ch4,alt,out.scaling_A_layer,out.geo.layer_dens.ch4,out.geo.air,out.geo.center_alts);
plot(ch4,alt,'r-','linewidth',2)
plot(ch4_smooth./out.geo.air',out.geo.center_alts,'k-','linewidth',2)
legend('Scaled prior','Aircore','Smoothed Aircore')
set(gca,'ylim',[0 30])

