
clear all
close all

% set path for absorption coeffs:
voigt_path = '/home/tukiains/Dropbox/voigt_shapes/';

[pathstr,name] = fileparts(which('get_ftir_files.m'));

% all files in input folder
%prefix = '/home/tukiains/Dropbox/ftir_spectra/';
prefix = [pathstr,'/../input_data/ftir_spectra/'];
%prefix = '/home/tukiains/Documents/ggg-2014/ggg-stable/i2s/opus-i2s/spectra/20130903/';

mfiles = dir([prefix,'so*']);
mfiles = struct2cell(mfiles);
mfiles = mfiles(1,:);

n = 5;
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

[co2,co2e,ch4,ch4e,co,coe,pres,alt,temp,air] = read_aircore_sounding(ac_file{1});
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

figure(4);
clf
ha = tight_subplot(2,1,[0.05],[0.15 0.01],[0.1 0.01])
axes(ha(1))
plot(out.wn,out.t,'linewidth',2,'color',[.3 .3 .3])
ylabel('FTIR spectrum')
set(ha(1),'xticklabel','')
% lower
axes(ha(2))
hold on
plot(out.wn,out.scaling_residual,'-','linewidth',2,'color', [0.2157 ...
                    0.4941 0.7216]);
plot(out.wn,out.dr_lm_residual(1:end-4),'-','linewidth',2,'color',[1.0000 ...
                      0.4980         0]);
set(ha,'xlim',[min(out.wn) max(out.wn)])
set(ha,'box','on')
grid on
xlabel('Wavenumber [cm^{-1}]')
ylabel('Residual')
set(ha(2),'ylim',[-4 2])
x = [6003:0.5:6005];
y = [-4:2:2];
set(ha,'xtick',x)
set(ha(2),'xticklabel',x)
set(ha(2),'ytick',y)
set(ha(2),'yticklabel',y)
h = legend('Prior scaling','Profile retrieval')
set(h,'location','southwest')
pos = get(h,'position');
pos(2) = pos(2)-0.01;
pos(1) = pos(1)+0.02;
set(h,'fontsize',11)
set(h,'position',pos)
legend boxoff
print_fig(15,10,'resis')