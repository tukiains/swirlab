function out = ch4_simu(voigt_root_path,mfile,k,per_alt,per)
% out = ch4_simu(voigt_root_path,mfile)
%

% date of the measurement
mdate = get_date(mfile);

% rootpath of swirlab
labpath = fileparts(which('calc_direct_geo.m'));

% select window 
[window,wnrange,gasvec,invgas,ninvgas,sol_shift_wn,solar_line_file,mindep] = window_details('ch4',3);

invgas = {'ch4','h2o'};
ninvgas = length(invgas);

% solar zenith angle
out.sza = get_sza_angle(mfile,[labpath,'/../input_data/ggg_runlog_files/so',mdate,'.grl'])

% voig path
voigtpath = [voigt_root_path,'voigt_shapes_',mdate,'_', ...
             int2str(min(window)),'-',int2str(max(window)),'/'];

% cross sections 
[c_wn,cros_o,c_alt] = read_cross_sections(gasvec, wnrange, labpath, voigtpath);

% ggg prior atmosphere file
afile = [labpath,'/../input_data/ggg_apriori_files/so',mdate,'.mav'];

% noise of the spectrum 
noise = 0.0015;

% reference (FTIR measurement)
[refe,wn] = read_ftir_spectrum(mfile,wnrange);

% altitude grid for the retrieval
alt = create_layering(70,100,1.002);

% direct sun geometry
[geo,cros] = calc_direct_geo(c_wn,cros_o,c_alt,wn,gasvec,afile,out.sza,alt);

% solar spectrum
sol = calc_solar_spectrum(wn,[labpath,'/../',solar_line_file]);
sol = sol./max(sol); % normalize (ftir spectrum not in physical units)

% lagrange matrix
wns = (wn-mean(wn))./std(wn);
wnref = wns([1,fix(end/2),end]);
L = lagrange(wnref,wns);

% default value for offset 
offset = 1e-4;

% remove some edges of the fitting window
ncut = 16;

% simulated measurement
ac_file = get_aircore_file(mdate,[labpath, '/../input_data/aircore/']);
[refe,wn_shift,sol_shift,simuatmos] = simulate_ftir_spectrum(c_wn,cros_o,c_alt,wn,gasvec,afile,out.sza,L,sol,noise,ac_file,per_alt,per);

% input for residual / jacobian calculation
varargin = create_varargin(wn,gasvec,cros,refe,invgas,sol,wn_shift,noise,L,geo,offset,ncut);

%% ----------------------------
%% prior scaling method with LM
%% ----------------------------

% prior profile of the first gas
x0 = geo.layer_dens.(char(invgas(1)));

% initial values
theta0 = [ones(1,ninvgas), 0.5, 0.5, 0.5, offset]; % retrieve also offset     

% LM-fit of scaling factors
[theta,cmat,rss,r] = levmar(@resfun,@jacfun,theta0,varargin);

% Averaging kernel (maybe incorrect way?)
nalt = length(geo.center_alts);
% use dimension reduction Jacobian function:
for n=1:ninvgas
    P(n) = {ones(nalt,1)};
end
theta2 = [log(theta(1:ninvgas)); theta(ninvgas+1:end)];
[~,~,K1] = jacfun_dr(theta2,ones(ninvgas,1),P,varargin);

[out.scaling_A_alpha,out.scaling_A_layer,out.scaling_A_column] = avek_scale(K1,x0,geo.los_lens,varargin{11},geo.altgrid,ncut);

% save results for output
out.geo = geo;
out.scaling_factors = theta;
out.scaling_residual = r;

% column
ind = find(geo.center_alts<20);
out.scaling_colu = trapz(geo.center_alts(ind),x0(ind));

%% ----------------------------------
%% prior covariance reduction with LM
%% ----------------------------------

% number of parameters
d = [k ones(1,ninvgas-1)];

% truncate prior covariance 
[P, C, Q] = reduce_dim(invgas,d,geo.center_alts);

% initial values
theta0 = [zeros(1,sum(d)) theta(ninvgas+1:end)'];
npar = length(theta0);

% jacobian and residual functions
jacfuni = @(theta0) jacfun_dr(theta0,d,P,varargin);
resfuni = @(theta0) resfun_dr(theta0,d,P,varargin);

% LM-fit of alpha parameters 
[theta2,cmat2,rss2,r2] = levmar(resfuni,jacfuni,theta0);

% averaging kernel
[~,~,J] = jacfuni(theta2);
[out.A_alpha,out.A_layer,out.A_column] = avek_dr(J,P{1},theta2(1:d(1)),x0,geo.los_lens,varargin{11},geo.altgrid,ncut);

% retrieved profiles etc. for output
out.dr_lm_atmos = redu2full(theta2,d,P,invgas,geo.layer_dens);
out.dr_lm_theta = theta2;
out.dr_lm_residual = r2;
out.dr_lm_cmat = cmat2;
out.dr_pri_C = C;
out.dr_lm_P = P;
out.dr_k = k;

% column
out.dr_colu = trapz(geo.center_alts(ind),out.dr_lm_atmos.ch4(ind));

