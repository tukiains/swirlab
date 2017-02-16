
clear all
close all

% set path for absorption coeffs:
voigt_path = '/home/tukiains/Dropbox/voigt_shapes/';

[pathstr,name] = fileparts(which('get_ftir_files.m'));

% all files in input folder
%prefix = [pathstr,'/../input_data/ftir_spectra/'];
prefix = '/media/tukiains/3b48ca75-37ff-42ba-ab71-b78f39cd9a79/tukiains/Documents/ggg-2014/ggg-stable/i2s/opus-i2s/spectra/20140508/';
mfiles = list_files(prefix,'so*');

% try some 
mfile = mfiles{fix(length(mfiles)/2)};

zenlim = 85;
usesimu = false;
lm_only = false;
lis = true;
k = 3;
% fixed offset?
fixo.scale = true; 
fixo.dr = true;
fixo.mcmc = true;

% retrieve ch4
out = ftir_dimred_mcmc(voigt_path,mfile,lm_only,lis,k,fixo,usesimu,zenlim);

air = out.geo.air;
alts = out.geo.center_alts;
x0 = out.geo.layer_dens.ch4;

% o2 retrieval is experimental
%out_o2 = o2_retrieval(voigt_path,mfile,zenlim);

figure(1)
clf
hold on
fa = 0.6;

% prior
if (lis)
    P = out.lis_P;
else
    P = out.dr_lm_P;
end
[h1 pmix] = show_prior(alts,x0,air,P{1});

set(h1,'facealpha',fa-0.2)

if (lm_only)
    % dimension reduction with LM:
    h2 = show_lm(alts,x0,air,out.dr_lm_P{1},out.dr_lm_theta,out.dr_lm_cmat,out.dr_k);
    set(h2,'facealpha',fa+0.2)
    str = 'DR (LM)';
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
%plot(out.geo.layer_dens.ch4./out.geo.air*out.scaling_factors(1),out.geo.center_alts,'b-','linewidth',2)

% aircore
ac_file = get_aircore_file(get_date(mfile),[pathstr, '/../input_data/aircore/']);

[co2,co2e,ch4,ch4e,co,coe,pres,alt,temp,air] = read_aircore_sounding(ac_file);
plot(ch4,alt,'r-','linewidth',2)

% smoothed aircore
[ch4_smooth aci] = smooth_ac(ch4,alt,out.A_layer,out.geo.layer_dens.ch4,out.geo.air,out.geo.center_alts);
plot(ch4_smooth,out.geo.center_alts,'k-','linewidth',2)

set(gca,'ylim',[0 70])
set(gca,'xlim',[0 2100])


