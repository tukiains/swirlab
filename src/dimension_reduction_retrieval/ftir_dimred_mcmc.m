
clear all
close all

% set path for the absorption coeffs:
% ----------------------------------------------------
voigt_root_path = '/home/tukiains/data/voigt_shapes/';
% ----------------------------------------------------

% select date yyyymmdd of the FTS measurement
% -------------------------------------------
mdate = '20140319';
%------------------

% rootpath of swirlab
labpath = fileparts(which('calc_direct_geo.m'));

% select measurement file (there might be several for this day)
mfiles = get_ftir_files(mdate,[labpath,'/../input_data/ftir_spectra/']);
mfile = cell2mat(mfiles(1)); % just pick one

% select window 
[window,wnrange,gasvec,invgas,sol_shift_wn,solar_line_file,mindep] = window_details('ch4',3);

% solar zenith angle
sza = get_sza_angle(mfile,[labpath,'/../input_data/ggg_runlog_files/so',mdate,'.grl']);

% voig path
voigtpath = [voigt_root_path,'voigt_shapes_',mdate,'_', ...
             int2str(min(window)),'-',int2str(max(window)),'/'];

% cross sections 
[c_wn, cros, c_alt] = read_cross_sections(gasvec, wnrange, labpath, voigtpath);

% ggg prior atmosphere file
afile = [labpath,'/../input_data/ggg_apriori_files/so',mdate,'.mav'];

% noise of the spectrum 
noise = 0.001;

% reference (FTIR measurement)
[refe,wn] = read_ftir_spectrum(mfile,wnrange);

% apodization
f = ggg_ils(2,length(wn),median(wn),median(diff(wn)),45,0);
f = f./sum(f);
refe = conv(refe,f,'same');

% altitude grid for the retrieval
alt = create_layering(70,70,1.01);

% direct sun geometry
[geo,cros] = calc_direct_geo(c_wn,cros,c_alt,wn,gasvec,afile,sza,alt);

% solar spectrum
sol = calc_solar_spectrum(wn,[labpath,'/../',solar_line_file]);
sol = sol./max(sol); % normalize (ftir spectrum not in physical units)

% lagrange matrix
wns = (wn-mean(wn))./std(wn);
wnref = wns([1,fix(end/2),end]);
L = lagrange(wnref,wns);

% wavelength shift
wn_shift = calc_wn_shift(geo,wn,gasvec,cros,refe,sol,mindep,L);

% solar wl shift
sol_shift = calc_sol_shift(wn,refe,sol,wn_shift,sol_shift_wn);    
sol = interp1(wn+sol_shift,sol,wn,'linear','extrap');

% input for residual / jacobian calculation
varargin = create_varargin(wn,gasvec,cros,refe,invgas,sol,wn_shift,noise,L,geo);

%% -------------------------
%% prior scaling method (LM)
%% -------------------------

% initial values
theta0 = [ones(1,length(invgas)), 0.5, 0.5, 0.5, 0.001];

% LM-fit of scaling factors
[theta,cmat,rss,r] = levmar(@resfun,@jacfun,theta0,varargin);

%% -------------------------------
%% dimension reduction method (LM)
%% -------------------------------

% number of parameters
d = [3 ones(1,length(invgas)-1)];

% truncate prior covariance 
[P, C, Q] = reduce_dim(invgas,d,geo.center_alts);

% initial values
theta0 = [zeros(1,sum(d)) theta(end-3:end)'];

% jacobian and residual functions
jacfuni = @(theta0) jacfun_dr(theta0,d,P,varargin);
resfuni = @(theta0) resfun_dr(theta0,d,P,varargin);

% LM-fit of alpha parameters 
[theta2,cmat2,rss2,r2] = levmar(resfuni,jacfuni,theta0);

% retrieved profiles
atmos2 = redu2full(theta2,d,P,invgas,geo.layer_dens);

%% ---------------------------------
%% dimension reduction method (MCMC)
%% ---------------------------------

simulations = 1000;                   % number of mcmc laps
updatesigma = 1;                      % 1 = estimate sigma2 based on residuals, 0 = constant / error assumed to be known
                                                                                                                                         
clear options params model

model.ssfun = @swirlab_ssfun_dimred;  % sum of squares function
model.sigma2 = 1;                     % initial error variance
model.N = length(refe);               % total number of observations
                                      % if update sigma is used
                                      % model.S20 = sigma2;                % prior for sigma2
                                      % model.N0 = 1;                      % prior accuracy for S20
options.savepostinss = 1;             % if 1 saves posterior ss, if 0 saves likelihood ss
options.nsimu = simulations;          % number of MCMC laps
options.burnintime = 5000;
options.waitbar = 1;                % graphical waitbar
options.verbosity = 5;              % how much to show output in Matlab window
options.printint = 100;             % how often to show info on acceptance ratios
options.updatesigma = updatesigma;  % 1 allow automatic sampling and estimation of the error variance
options.adaptint = 100;              % how often to adapt the proposal
options.method = 'am';              % adaptation method: 'mh', 'dr', 'am' or 'dram'
options.stats = 1;
options.qcov = eye(npar)*0.1^2;     % initial covariance for the Gaussian proposal density of the MCMC sampler
                                    % options.qcov = diag(diag(cmat));

% MCMC parameters to sample: lognormal prior with mode 1, variance 1
% (sqrt(0.32228)), mean 1.62161 (0.32228) mu = data.mu;

% alpha-parameters
for i=1:sum(d)
    %            'name',              initial value,    min,    max,    prior mean,     prior std deviation
    %params{i} = {num2str(i),          theta2(i),      -Inf, Inf,    0,              1};
    params{i} = {num2str(i),          0,      -Inf, Inf,    0,              1};
end

% parametrized polynomial terms
for i=sum(d)+1:npar
    params{i} = {num2str(i),          theta2(i),        0,   2,    0.2,             1};
end

% offset term
params{end} = {num2str(npar),     0.001,   0,  1,    0.01,           Inf};

data.type = 'mcmc';

% Run MCMCM
[res,chain,s2chain,sschain] = mcmcrun(model,data,params,options);


% show result
close all
[P, C, Q] = reduce_dim(invgas, gasvec, d, alt2, layer_dens, airi, samples);
figure(1)
hold on
gas_prior = interp1(alt,data.atmos(1,:)./air,alt2);
skip = 1;
redchain = chain(end-20000:end,:);
Proj = P{1};
A = redchain(:,1:d(1))*Proj';

profchain = exp(bsxfun(@plus,A,log(gas_prior.*airi))); %log-normal case
%profchain = bsxfun(@plus,A,gas_prior.*airi); % linear case
profs = bsxfun(@rdivide,profchain,airi);
%profchain = bsxfun(@plus,A,gas_prior.*airi);
%profs = bsxfun(@rdivide,profchain,airi);

plot(atmos(1,:)./air*1e9,alt,'b--','linewidth',2)
predplot_simo(alt2,plims(profs*1e9,[0.025 0.5 0.975]),[.5 .7 .3],[.5 .7 .3]);
plot(ac_ch4,ac_alt,'r-','linewidth',2)
%plot(ac_ch4,hpa2km(ac_pres),'m-','linewidth',2)

h = legend('Draws from prior','Prior mean','95% posterior envelope','AirCore')
legend boxoff
set(gca,'xlim',[0.4e-6 2.3e-6]*1e9)
ylabel('Altitude [km]','fontsize',14)
set(gca,'ylim',[0 40])
set(gca,'box','on')
title('Sodankyla, Finland, 3.9.2013','fontsize',14)
set(gca,'xscale','lin')
set(h,'location','southwest')
set(h,'fontsize',14)
set(gca,'xtick',[400:300:2200])
xlabel('CH_4 mixing ratio [ppb]','fontsize',14)
xlabh = get(gca,'XLabel');
set(xlabh,'Position',get(xlabh,'Position')-[0 .2 0])
%print_fig(18,16,'mcmc_dimred_3_9_2013')

if (1==2)
    figure(2)
    clf
    res.names = {'ch4-1','ch4-2','ch4-3','p-1','p-2','p-3','offset'}
    mcmcplot(redchain,[],res,'chainpanel');
    print_fig(22,18,'mcmc_dimred_chain')
    
    figure(3)
    %res.names = {'','','','','','',''}
    clf
    mcmcplot(redchain(end-20000:end,:),[],res,'pairs');
    nvec = [1,7,13,19,25,31]
    for n=nvec
        subplot(6,6,n)
        set(gca,'ytick',[])
        set(gca,'yticklabel',{})
    end
    nvec = [31:36]
    for n=nvec
        subplot(6,6,n)
        set(gca,'xtick',[])
        set(gca,'xticklabel',{})
    end
    print_fig(18,18,'mcmc_dimred_pairs')
end

%figure(3)
%mcmcplot(redchain,[],res,'pairs', 2);

% res has many details about the mcmc run
% chain is the actual MCMC chain for the parameters
% s2chain is the chain for the observation error variance
% sschain has the chain of the ssfun results

%% SAVES DATA
RES.res = res;
RES.chain = chain;
RES.s2chain = s2chain;
RES.sschain = sschain;
%save('/home/tukiains/Documents/MATLAB/swirlab/dimension_reduction_retrieval/mcmc/results/3-9-2013-100k_10kmcorr.mat','RES','alt','data','alt2','P','d','gas_prior','airi','ac_ch4','ac_alt','air')

% column estimate
colu = posterior_column(profs(end-10000:end,:),airi,alt2,20)