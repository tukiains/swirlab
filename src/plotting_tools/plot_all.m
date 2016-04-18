
% requires tight_subplot.m
% available at: http://goo.gl/0BD7wA

clear all
close all

data = load('/home/tukiains/Dropbox/Public/ftir_results.mat');
%data = load('/home/tukiains/Dropbox/Public/ftir_results_lm.mat');
%data = load('/home/tukiains/Dropbox/Public/ftir_results_lis_k=6.mat');

ac_path = '/home/tukiains/Documents/MATLAB/swirlab_git/swirlab/input_data/aircore/';

ha = tight_subplot(5,2,[.02 .01],[.15 .01],[.1 .01])

scale = 1e9;

for n = 1:10;
    
    figure(1)
    out = data.out(n);
    axes(ha(n));

    %h2 = show_lm(out.geo.center_alts,out.geo.layer_dens.ch4*scale,out.geo.air,out.dr_lm_P{1},out.dr_lm_theta,out.dr_lm_cmat,out.dr_k);
    h2 = plot_curtain(out.geo.center_alts,plims(out.mcmc_profs*scale,[0.025 0.5 0.975]),[.5 .7 .3]);    

    set(h2,'facealpha',1)
   
    hold on
   
    plot(out.geo.layer_dens.ch4./out.geo.air*1e9,out.geo.center_alts,'b--','linewidth',2)
    
    dd = get_date(data.mfile(n,:));

    % aircore
    ac_file = get_aircore_file(dd,ac_path);
    [co2,co2e,ch4,ch4e,co,coe,pres,alt,temp,air] = read_aircore_sounding(ac_file);
    plot(ch4,alt,'r-','linewidth',2) 

    %smoothed version
    [ch4_smooth, aci] = smooth_ac(ch4/1e9,alt,out.A_layer,out.geo.layer_dens.ch4,out.geo.air,out.geo.center_alts,1);
    plot(ch4_smooth*1e9./out.geo.air',out.geo.center_alts,'k-','linewidth',2)

    ddd = datestr(datenum([dd(7:8),'.',dd(5:6),'.',dd(1:4)],'dd.mm.yyyy'),1);

    text(300,5,ddd,'fontsize',10)
    
    if (n>=9)        
        xlabel('CH_4 [ppb]','fontsize',10)        
    end

    if (mod(n,2)==1)
        ylabel('Altitude [km]','fontsize',10)
    end

end

set(ha,'ylim',[0 30])
set(ha,'xlim',[200 2100])

x = 400:400:2000;
y = 0:10:30;

set(ha,'xtick',x);
set(ha,'ytick',y);
set(ha,'yticklabel','');
set(ha,'xticklabel','');

set(ha(9),'xticklabel',x,'fontsize',10);
set(ha(10),'xticklabel',x,'fontsize',10);

for n=1:2:9
    set(ha(n),'yticklabel',y,'fontsize',10);
end

l = legend('95% Posterior','Prior mean','AirCore','AirCore (smoothed)')
set(l,'position',[0.5 0.05 0.1 0.02])
set(l,'fontsize',10)
legend boxoff

%print_fig(15,21,'mcmc_profiles')


