
% requires tight_subplot.m
% available at: http://goo.gl/0BD7wA

clear all
close all

data = load('/home/tukiains/Dropbox/Public/ftir_results.mat');

ac_path = '/home/tukiains/Documents/MATLAB/swirlab/input_data/aircore/';

ha = tight_subplot(5,2,[.01 .03],[.01 .01],[.01 .01])

for n = 1:10;
    
    out = data.out(n);
    axes(ha(n));

    %h2 = show_lm(out.geo.center_alts,out.geo.layer_dens.ch4,out.geo.air,out.dr_lm_P{1},out.dr_lm_theta,out.dr_lm_cmat,out.dr_k);

    h2 = plot_curtain(out.geo.center_alts,plims(out.mcmc_profs,[0.025 0.5 0.975]),[.5 .7 .3]);

    set(h2,'facealpha',1)
   
    hold on
    
    ac_file = get_aircore_file(get_date(data.mfile(n,:)),ac_path);
    [co2,co2e,ch4,ch4e,co,coe,pres,alt,temp,air] = read_aircore_sounding(ac_file);
    plot(ch4./1e9,alt,'r-','linewidth',2)
    
end

set(ha,'ylim',[0 30])
set(ha,'xlim',[0.2e-6 2.1e-6])
set(ha,'yticklabel','');
set(ha,'xticklabel','')

%print_fig(10,18,'mcmc_profiles')


