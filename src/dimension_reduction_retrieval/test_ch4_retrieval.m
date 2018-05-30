clear all
close all

% set SWIRLAB root path
addpath(genpath('/home/tukiains/Temp/swirlab'))

% set path to Marko Laine's MCMC Toolbox
addpath(genpath('/home/tukiains/Documents/MATLAB/mcmcstat'))

% set path for absorption coeffs:
voigt_path = '/home/tukiains/Dropbox/voigt_shapes/';

[pathstr,name] = fileparts(which('get_ftir_files.m'));

% all files in input folder
prefix = [pathstr,'/../input_data/ftir_spectra/'];

mfiles = list_files(prefix,'so*');

% options
zenlim = 82;        % limit for solar zenith angle
usesimu = false;    % use simulated measurement
lm_only = false;    % optimal estimation retrieval without MCMC
lis = true;         % use LIS
ace = true;         % use empirical ACE-prior
jaco_sample = true; % evaluate LIS basis at prior mean if false, sample from levmar-approximation if true
k = 3;              % dimension of prior reduction and LIS subspace, optimal k = 3

% fixed offset?
fixo.scale = true; 
fixo.dr = true;
fixo.mcmc = true;

% select file
mfile = mfiles{2};

% retrieve ch4
out = ftir_dimred_mcmc(voigt_path,mfile,lm_only,lis,k,fixo,usesimu,zenlim,jaco_sample,ace);

% physical parameters for plotting
air = out.geo.air;
alts = out.geo.center_alts;
x0 = out.geo.layer_dens.ch4;
wn = out.wn;

% o2 retrieval is experimental
%out_o2 = o2_retrieval(voigt_path,mfile,zenlim);

figure(1)
clf
hold on
fa = 0.6;

% prior
P = out.dr_lm_P;

if lis
    P2 = out.full_P;
    [h1 pmix] = show_prior(alts,x0,air,P2{1});
else 
    [h1 pmix] = show_prior(alts,x0,air,P{1});
end 

set(h1,'facealpha',fa-0.2)

if (lm_only)
    % dimension reduction with LM:
    h2 = show_lm(alts,x0,air,P{1},out.dr_lm_theta,out.dr_lm_cmat,out.dr_k);
    set(h2,'facealpha',fa+0.2)
    
    if (lis)
        str = 'LIS (LM)';
    else
        str = 'DR (LM)';
    end
else
    % dimension reduction with MCMC (2-sigma posterior):
    h2 = plot_curtain(alts,plims(out.mcmc_profs,[0.025 0.5 0.975]),[.5 .7 .3]);
    set(h2,'facealpha',fa+0.2)
    if (lis)
        str = 'LIS (MCMC)';
    else
        str = 'DR (MCMC)';
    end
end

% prior
plot(x0./air*1e9,alts,'b--','linewidth',2)

% scaled prior
%plot(out.geo.layer_dens.ch4./out.geo.air*out.scaling_factors(1)*1e9,out.geo.center_alts,'k-','linewidth',2)

% aircore
ac_file = get_aircore_file(get_date(mfile),[pathstr, '/../input_data/aircore/']);

[co2,co2e,ch4,ch4e,co,coe,pres,alt,temp,air] = read_aircore_sounding(ac_file);
plot(ch4,alt,'r-','linewidth',2)

% smoothed aircore
[ch4_smooth aci] = smooth_ac(ch4,alt,out.A_layer,out.geo.layer_dens.ch4,out.geo.air,out.geo.center_alts);
plot(ch4_smooth,out.geo.center_alts,'k-','linewidth',2)

set(gca,'ylim',[0 70])
set(gca,'xlim',[0 2100])
ylabel('Altitude [km]')
xlabel('CH4 [ppb]')
title([str, ', k: ', num2str(k)]);

% model fit and residual plot
figure(2)
clf
subplot(2,1,1)
hold on
plot(wn,out.t,'k-')
plot(wn,out.dr_fitted_model,'ro','markerfacecolor','r','markersize',4)
set(gca,'xlim',[min(wn) max(wn)])
grid on
ylabel('Transmittance')
xlabel('Wavenumber [cm^{-1}]')
h = legend('Measurement','Model at optimum');
set(h,'location','southeast')
legend boxoff
subplot(2,1,2)
hold on
plot(wn,out.dr_lm_residual(1:length(wn)),'k-')
set(gca,'xlim',[min(wn) max(wn)])
grid on
ylabel('Residual')
xlabel('Wavenumber [cm^{-1}]')