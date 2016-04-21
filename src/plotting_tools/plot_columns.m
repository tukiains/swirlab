
clear all
close all

%mdate = '20130903';
%mdate = '20130905';
mdate = '20131022';
%mdate = '20140319';
%mdate = '20140409';
%mdate = '20140508';
%mdate = '20140714';
%mdate = '20140715';
%mdate = '20140716';
%mdate = '20140902';
%mdate = '20141105'; 


% read swirlab results
data = load(['/home/tukiains/data/dimension_reduction_results/retrieved_offset/ftir_results_',mdate,'.mat']);

col_file = ['/home/tukiains/data/opus/ggg_results/',mdate,'/ch4_6002.so',mdate,'.col'];
ak_path = ['/home/tukiains/data/opus/ggg_results/',mdate,'/aks/ch4_6002/'];
grl_file = ['/home/tukiains/data/opus/ftir_igrams/aircore/grl_files/so',mdate,'.grl'];

% read ggg columns etc
[ggg_col_name, ggg_sza, ggg_VSF_ch4, ggg_OVC_ch4, ggg_VSF_h2o, ggg_OVC_h2o] = read_col(col_file,grl_file);

% read ggg AKs
[ggg_ak_name, ggg_ak, ggg_ak_level] = read_ggg_avek(ak_path);
ggg_ak_level = double(ggg_ak_level);

% our zenith angles
sza = cat(1,data.out.sza);

del = 1e-12;

for n=1:length(sza)

    sf(n) = data.out(n).scaling_factors(1);

    ak(n,:) = data.out(n).scaling_A_column./diff(data.out(1).geo.altgrid);

    ggg_ind(n) = find(ggg_sza==sza(n),1);

end

% just take ggg data that we have also
ggg_ak_name = ggg_ak_name(ggg_ind,:);
ggg_ak = ggg_ak(ggg_ind,:);
ggg_VSF_ch4 = ggg_VSF_ch4(ggg_ind);
ggg_OVC_ch4 = ggg_OVC_ch4(ggg_ind);
ggg_sza = ggg_sza(ggg_ind);

alt = data.out(1).geo.center_alts; % altitude vector
x0 = data.out(1).geo.layer_dens.ch4; % prior
air = data.out(1).geo.air; % air density

alt_dense = [0:0.01:70];
x0i = exp(interp1(alt,log(x0),alt_dense,'linear','extrap'));
airi = exp(interp1(alt,log(air),alt_dense,'linear','extrap'));
ovc_simo = trapz(alt_dense,x0i)*1e5;


figure(1)
clf
hold on
plot(sza,sf,'ro')
plot(ggg_sza,ggg_VSF_ch4,'bo');
set(gca,'xlim',[min(sza) max(sza)])
legend('SWIRLAB scaling factor','GGG scaling factor')
legend boxoff
%print_fig(15,10,'scaling_factor_comparison')

figure(2);
clf
subplot(1,2,1)
plot(ak,alt,'r-')
set(gca,'xlim',[0 1.5])
legend('SWIRLAB')
legend boxoff
subplot(1,2,2)
plot(ggg_ak,ggg_ak_level,'b-')
set(gca,'xlim',[0 1.5])
legend('GGG')
legend boxoff
%print_fig(15,10,'ak_comparison')

%% calculate columns
ac_file = get_aircore_file(mdate,'/home/tukiains/Documents/MATLAB/swirlab_git/swirlab/input_data/aircore/');
[~,~,ac_ch4,~,~,~,~,ac_alt,~,~] = read_aircore_sounding(ac_file);
ac_ch4 = ac_ch4/1e9;

% extrapolate aircore
aci = extrapolate_ac(ac_ch4,ac_alt,airi,alt_dense);
aci2 = extrapolate_ac(ac_ch4,ac_alt,air,alt);

ggg_aki = interp1(ggg_ak_level,ggg_ak',alt_dense);

swirlab_aki = interp1(alt,ak',alt_dense,'linear','extrap');

ind = find(alt_dense>0 & alt_dense<20);

for n=1:length(sza)

    % ac col with ggg AK
    p1 = aci.*ggg_aki(:,n)';
    ac_col_ggg(n) = trapz(alt_dense(ind),p1(ind))*1e5; % this
    
    % ac col with swirlab AK
    p1 = aci.*swirlab_aki(:,n)';
    ac_col_swirlab(n) = trapz(alt_dense(ind),p1(ind))*1e5; % this

    % This is another way. gives ~ the same result as above.
    %akf = data.out(n).scaling_A_layer;
    %ac_smooth = x0' + akf*(aci2-x0)';
    %ac_smooth = exp(interp1(alt,log(ac_smooth),alt_dense,'linear','extrap'));
    %ac_col_swirlab2(n) = trapz(alt_dense, ac_smooth)*1e5;

    % ggg col. these two are approximately the same
    p1 = x0i*ggg_VSF_ch4(n);
    ggg_col(n) = trapz(alt_dense(ind),p1(ind))*1e5; % this
    %ggg_col2(n) = ggg_VSF_ch4(n)*ggg_OVC_ch4(n);

    % swirlab col (with scaling)
    p1 = x0i*sf(n);    
    swirlab_col(n) = trapz(alt_dense(ind),p1(ind))*1e5; % this

    % dimension reduction solution....
    akf = data.out(n).A_layer;
    ac_smooth = x0' + akf*(aci2-x0)';
    ac_smooth = exp(interp1(alt,log(ac_smooth),alt_dense,'linear','extrap'));

    ac_col_swirlab2(n) = trapz(alt_dense(ind), ac_smooth(ind))*1e5; % this

    prof = data.out(n).dr_lm_atmos.ch4;
    profi = exp(interp1(alt,log(prof),alt_dense,'linear','extrap'))*1e5;    
 
    dr_col(n) = trapz(alt_dense(ind),profi(ind)); % this


end

figure(3)
clf
hold on

% swirlab scaling
plot(sza,ac_col_swirlab,'b-o','markerfacecolor','b')
plot(sza,swirlab_col,'b-*')

%plot(sza,ac_col_swirlab2,'m*')

% ggg scaling
%plot(sza,ac_col_ggg,'go')
%plot(sza,ggg_col,'g*')

%plot(sza,ggg_col2,'r*')

% swirlab DR
plot(sza,ac_col_swirlab2,'r-o','markerfacecolor','r')
plot(sza,dr_col,'r-*')

set(gca,'xlim',[min(sza) max(sza)])
%h = legend('AirCore (swirlab scaling AK)','swirlab (scaling)','AirCore (GGG AK)','GGG','AirCore (DR AK)','DR-LM')
legend boxoff
set(h,'location','northwest')

title(mdate)
ylabel('CH4 column (0-20km)')
xlabel('SZA')
set(gca,'box','on')

print_fig(22,15,'cols_fixo')

