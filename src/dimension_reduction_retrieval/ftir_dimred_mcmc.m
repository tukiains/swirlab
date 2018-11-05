function out = ftir_dimred_mcmc(voigt_root_path,mfile,lm_only,lis,k,fixoff,usesimu,zenlim,jaco_sample,ace)
% out = ftir_dimred_mcmc(voigt_root_path,mfile,lm_only,lis,k,fixoff,usesimu,zenlim,jaco_sample,ace)
%

% date of the measurement
mdate = get_date(mfile);

% rootpath of swirlab
labpath = fileparts(which('calc_direct_geo.m'));

% select window 
[window,wnrange,gasvec,invgas,ninvgas,sol_shift_wn,solar_line_file] = window_details('ch4',3);

% solar zenith angle
out.sza = get_sza_angle(mfile,[labpath,'/../input_data/ggg_runlog_files/so',mdate,'.grl']);

if (isempty(out.sza) | out.sza>zenlim)
    disp('Error: cant find solar zenith angle for this scan (or too high)')
    return
end

voigtpath = [voigt_root_path,'voigt_shapes_',mdate,'_', ...
             int2str(min(window)),'-',int2str(max(window)),'/'];

% cross sections 
[c_wn,cros_o,c_alt] = read_cross_sections(gasvec, wnrange, labpath, voigtpath);

% ggg prior atmosphere file
afile = [labpath,'/../input_data/ggg_apriori_files/so',mdate,'.mav'];

% noise of the spectrum 
init_noise = 0.0015;
noise = init_noise;

% reference (FTIR measurement)
[refe,wn] = read_ftir_spectrum(mfile,wnrange);

if (isempty(refe))
    msg = 'Error: No spectrum';
    disp('Error: No spectrum')
    error(msg)
end

% apodization
f = ggg_ils(2,length(wn),median(wn),median(diff(wn)),45,0);
refe = conv(refe,f,'same');

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

% remove some edges of the fitting window
ncut = 16;

% wavelength shift
wn_shift = calc_wn_shift(geo,wn,gasvec,cros,refe,sol,L,ncut);

% solar wl shift
sol_shift = calc_sol_shift(mfile,wn_shift,labpath,solar_line_file,sol_shift_wn);
sol = interp1(wn+sol_shift,sol,wn,'linear','extrap');

% default value for offset 
offset = 0;
if (usesimu) % then simulate measurement
    ac_file = get_aircore_file(mdate,[labpath, '/../input_data/aircore/']);
    [refe,wn_shift,sol_shift] = simulate_ftir_spectrum(c_wn,cros_o,c_alt,wn,gasvec,afile,out.sza,L,sol,noise,ac_file);
end

% input for residual / jacobian calculation
varargin = create_varargin(wn,gasvec,cros,refe,invgas,sol,wn_shift,noise,L,geo,offset,ncut,false);

%% ----------------------------
%% prior scaling method with LM
%% ----------------------------

% prior profile of the first gas
x0 = geo.layer_dens.(char(invgas(1)));

% initial values
if (fixoff.scale)
    theta0 = [ones(1,ninvgas), 0.5, 0.5, 0.5]; % fixed offset
else
    theta0 = [ones(1,ninvgas), 0.5, 0.5, 0.5, offset]; % retrieve also offset     
end

% LM-fit of scaling factors
[theta,cmat,rss,r] = levmar(@resfun,@jacfun,theta0,varargin);

% Averaging kernel. Use dimension reduction Jacobian function.
nalt = length(geo.center_alts);
for n=1:ninvgas
    P(n) = {ones(nalt,1)};
end
theta2 = [log(theta(1:ninvgas)); theta(ninvgas+1:end)];
[~,~,K1] = jacfun_dr(theta2,ones(ninvgas,1),P,0,varargin);
[out.scaling_A_alpha,out.scaling_A_layer,out.scaling_A_column] = avek_scale(K1,x0,geo.los_lens,varargin{11},geo.altgrid,ncut);

% save results for output
out.geo = geo;
out.scaling_factors = theta;
out.scaling_residual = r;
out.wn = wn(ncut:end-ncut);
out.t = refe(ncut:end-ncut);
out.sol = sol(ncut:end-ncut);

nobs = length(r); % number of true wavelengths

%% ----------------------------------
%% prior covariance reduction with LM
%% ----------------------------------

% number of parameters
d = [k ones(1,ninvgas-1)];

% truncate prior covariance 
[P, C, Q] = reduce_dim(invgas,d,geo.center_alts, ace);

% initial values
theta0 = [zeros(1,sum(d)) theta(ninvgas+1:end)'];
if (fixoff.scale & not(fixoff.dr))
    theta0 = [theta0 offset];
elseif (not(fixoff.scale) & fixoff.dr)
    theta0 = theta0(1:end-1);
end

npar = length(theta0);

% jacobian and residual functions
O = 0;
jacfuni = @(theta0) jacfun_dr(theta0,d,P,O,varargin);
resfuni = @(theta0) resfun_dr(theta0,d,P,O,varargin);

% LM-fit of alpha parameters 
[theta2,cmat2,rss2,r2] = levmar(resfuni,jacfuni,theta0);

% update error estimate and retrieve again 
noise = noise*sqrt(rss2);
varargin = create_varargin(wn,gasvec,cros,refe,invgas,sol,wn_shift,noise,L,geo,offset,ncut,false);
jacfuni = @(theta0) jacfun_dr(theta0,d,P,O,varargin);
resfuni = @(theta0) resfun_dr(theta0,d,P,O,varargin);
[theta2,cmat2,rss2,r2] = levmar(resfuni,jacfuni,theta0);

% averaging kernel
[~,~,J] = jacfuni(theta2);
[out.A_alpha,out.A_layer,out.A_column] = avek_dr(J,P{1},theta2(1:d(1)),x0,geo.los_lens,varargin{11},geo.altgrid,ncut,geo.air);

% retrieved profiles etc. for output
out.dr_lm_atmos = redu2full(theta2,d,P,O,invgas,geo.layer_dens,geo.air,false);
out.dr_lm_theta = theta2;
out.dr_lm_residual = r2;
out.dr_lm_cmat = cmat2;
out.dr_pri_C = C;
out.dr_lm_P = P;
out.dr_k = k;

% model at optimum
err = cell2mat(varargin(11));
out.dr_fitted_model = out.t + r2(1:end-sum(d)).*err(ncut:end-ncut);
out.version = swirlab_version();
out.date = mdate;

%% --------------------------
%%       LIS method
%% --------------------------

if (lis)
    disp('using likelihood-informed dimension reduction')
    
    % Get the full-dimensional prior mean for jacobian
    d2 = [100 ones(1,ninvgas-1)]; 
    [P2, C2, Q2] = reduce_dim(invgas,d2,geo.center_alts, ace); 
    theta00 = [zeros(1,sum(d2)) theta(ninvgas+1:end)']; 
    if (fixoff.scale & not(fixoff.dr))
        theta00 = [theta00 offset];     
    elseif (not(fixoff.scale) & fixoff.dr)
        theta00 = theta00(1:end-1); 
    end
    npar2 = length(theta00); 
    noise = init_noise;
    varargin = create_varargin(wn,gasvec,cros,refe,invgas,sol,wn_shift,noise,L,geo,offset,ncut,false);
    jacfunii = @(theta00) jacfun_dr(theta00,d2,P2,O,varargin);
    
    out.full_P = P2; % for plotting prior in LIS case, fixme
    
    % measurement error
    Ly = diag(ones(nobs,1)./noise); % Cy = inv(Ly'*Ly)
    
    % covariance
    cov1 = C2{1};

    meanupd = @(x,m,i) m + 1./i*(x-m);
    Jm = 0;
        
    % cholesky of prior covariance
    Lx = inv(chol(cov1+eye(size(cov1))*0.1,'lower')); % Cx = inv(Lx'*Lx);
    if jaco_sample % sample LIS basis from LevMar -approximation
      disp('Sampling...')
      for i=1:500 % sufficient amount of samples needs to be investigated
        [~,~,Jsample] = jacfuni(mvnorr(1,theta2,cmat2)');
        Jsample = Jsample*diag(geo.air/1e9);
        Js = (Ly*Jsample/Lx)'*(Ly*Jsample/Lx);
        Jm = meanupd(Js,Jm,i);
      end
      disp('...done')
    else % evaluate LIS basis in prior mean
      [~,~,Jsample] = jacfunii(theta00');
      Jsample = Jsample*diag(geo.air/1e9);
      Jm = (Ly*Jsample/Lx)'*(Ly*Jsample/Lx);
    end
    
    % svd of that
    [~,s,v] = svd(Jm,0);
    
    % how many components is needed?
    %     - this would show the optimal k from svd:
    % k = length(find(diag(s)>1))+1;
    % d(1) = k;
    % npar = npar - (k_old-k);
    
    % projection matrix
    P(1) = {Lx\v(:,1:k)};
    
    % complement space
    Po = Lx\v(:,k+1:end);

    % projection from full to LIS space
    O = (Lx'*v(:,1:k))';
    
    % projection from full to CS
    Oo = (Lx'*v(:,k+1:end))';
    
    out.Lx = Lx;
    out.lis_s = diag(s);
    out.lis_P = P(1);
    varargin{14} = true;
    
    % --------------------- LIS OE --------------------- %
    
    theta0 = [zeros(1,sum(d)) theta(ninvgas+1:end)'];
    if (fixoff.scale & not(fixoff.dr))
        theta0 = [theta0 offset];
    elseif (not(fixoff.scale) & fixoff.dr)
        theta0 = theta0(1:end-1);
    end
    
    npar = length(theta0);
    
    % jacobian and residual functions
    jacfuni = @(theta0) jacfun_dr(theta0,d,P,O,varargin);
    resfuni = @(theta0) resfun_dr(theta0,d,P,O,varargin);
    
    % LM-fit of alpha parameters 
    [theta2,cmat2,rss2,r2] = levmar(resfuni,jacfuni,theta0);

    % update error estimate and retrieve again 
    noise = noise*sqrt(rss2);
    varargin = create_varargin(wn,gasvec,cros,refe,invgas,sol,wn_shift,noise,L,geo,offset,ncut,false);
    jacfuni = @(theta0) jacfun_dr(theta0,d,P,O,varargin);
    resfuni = @(theta0) resfun_dr(theta0,d,P,O,varargin);
    [theta2,cmat2,rss2,r2] = levmar(resfuni,jacfuni,theta0);

    % averaging kernel
    [~,~,J] = jacfuni(theta2);
    [out.A_alpha,out.A_layer,out.A_column] = avek_dr(J,P{1},theta2(1:d(1)),x0,geo.los_lens,varargin{11},geo.altgrid,ncut,geo.air);

    % retrieved profiles etc. for output
    out.dr_lm_atmos = redu2full(theta2,d,P,O,invgas,geo.layer_dens,geo.air,lis);
    out.dr_lm_atmos.ch4 = out.dr_lm_atmos.ch4+(Po*Oo*(x0'./geo.air'*1e9))'.*geo.air/1e9;
    out.dr_lm_theta = theta2;
    out.dr_lm_residual = r2;
    out.dr_lm_cmat = cmat2;
    out.dr_pri_C = C;
    out.dr_lm_P = P;
    out.dr_k = k;

    % model at optimum
    err = cell2mat(varargin(11));
    out.dr_fitted_model = out.t + r2(1:end-sum(d)).*err(ncut:end-ncut);
    out.version = swirlab_version();
    
else
    disp('using ordinary dimension reduction')
end

if (lm_only)
    return
end

%% -------------
%% MCMC sampling 
%% -------------

model.ssfun = @ssfun_mcmc;            % sum of squares function
model.sigma2 = 1;                     % initial error variance
model.N = nobs;                       % number of observations (true wavelengths)
options.savepostinss = 1;             % if 1 saves posterior ss, if 0 saves likelihood ss
options.nsimu = 40000;                % number of MCMC laps
options.burnintime = 1500;
options.waitbar = 1;                  % graphical waitbar
options.verbosity = 5;                % how much to show output in Matlab window
options.printint = 100;               % how often to show info on acceptance ratios
options.updatesigma = 1;              % 1 allow automatic sampling and estimation of the error variance
options.adaptint = 100;               % how often to adapt the proposal
options.method = 'am';                % adaptation method: 'mh', 'dr', 'am' or 'dram'
options.stats = 1;
options.qcov = diag(diag(C{1}));     % initial covariance for the Gaussian proposal density of the MCMC sampler

% MCMC parameters to sample: 

% alpha-parameters
for i = 1:sum(d)
    % 'name', initial value, min, max, prior mean, prior std
    params{i} = {num2str(i), 0, -Inf, Inf, 0, 1};
end

% parametrized polynomial terms
g = 2;
for i = sum(d)+1:npar
    params{i} = {num2str(i), theta2(end-g), 0, 1, 0.2, 1};
    g = g-1;
end

% offset term
if (fixoff.dr & not(fixoff.mcmc)) 
    npar = npar+1;
    params{end+1} = {num2str(npar), 0.001, 0, 1, 0.001, 0.5};
    options.qcov = diag([diag(cmat2);1e-4]);
elseif (not(fixoff.dr) & fixoff.mcmc) 
    npar = npar-1;
    params = params(1:end-1);
    options.qcov = diag(diag(cmat2(1:end-1,1:end-1)));    
elseif (not(fixoff.dr) & not(fixoff.mcmc))
    params{end} = {num2str(npar), 0.001, 0, 1, 0.001, 0.5};
end

data.d = d;
data.P = P;
data.O = O;
data.varargin = varargin;
data.x0 = x0;

% run MCMC
[res,chain,s2chain,sschain] = mcmcrun(model,data,params,options);

% retrieved profiles
redchain = chain(end-20000:end,:); % take last 20k samples?

fullchain = redchain(:,1:d(1))*P{1}' + x0(:)'./geo.air*1e9;

if (lis)
    fullchain = fullchain + randn(size(fullchain,1),size(Po,2))*Po';
end

% save for output
out.mcmc_profs = fullchain;
out.mcmc_res = res;
out.mcmc_chain = chain;
out.mcmc_s2chain = s2chain;
out.mcmc_sschain = sschain;
out.data = data;
