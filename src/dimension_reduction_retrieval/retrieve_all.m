function out = retrieve_all(voigt_root_path,mfiles,k,fixoff,zenlim)
% retrieve_all(voigt_root_path,mfiles,k,fixoff,zenlim)
%

% date of the measurement (they all should be from the same day!)
mdate = get_date(cell2mat(mfiles(1)));

% rootpath of swirlab
labpath = fileparts(which('calc_direct_geo.m'));

% select window 
[window,wnrange,gasvec,invgas,ninvgas,sol_shift_wn,solar_line_file,mindep] = window_details('ch4',3);

% voigt path
voigtpath = [voigt_root_path,'voigt_shapes_',mdate,'_', ...
             int2str(min(window)),'-',int2str(max(window)),'/'];

% cross sections 
[c_wn,cros_o,c_alt] = read_cross_sections(gasvec, wnrange, labpath, voigtpath);

% ggg prior atmosphere file
afile = [labpath,'/../input_data/ggg_apriori_files/so',mdate,'.mav'];

% retrieve all
jj = 1;
for ii=1:length(mfiles)

    mfile = cell2mat(mfiles(ii));

    % solar zenith angle
    sza = get_sza_angle(mfile,[labpath,'/../input_data/ggg_runlog_files/so',mdate,'.grl']);

    if (isempty(sza) | sza>zenlim)
        disp('Error: cant find solar zenith angle for this scan (or too high)')
        continue
    end

    % initial noise of the spectrum 
    noise = 0.0015;
    
    % reference (FTIR measurement)
    disp('Reading spectrum..')
    [refe,wn] = read_ftir_spectrum(mfile,wnrange);
    
    if (isempty(refe))
        disp('Error: No spectrum')
        continue
    end
    
    % apodization
    clear conv_spectrum
    f = ggg_ils(2,length(wn),median(wn),median(diff(wn)),45,0);
    refe = conv(refe,f,'same');
    
    % altitude grid for the retrieval
    alt = create_layering(70,100,1.002);
    
    disp('Calculating geometry..')
    % direct sun geometry
    [geo,cros] = calc_direct_geo(c_wn,cros_o,c_alt,wn,gasvec,afile,sza,alt);
    
    disp('Reading solar spectrum..')
    % solar spectrum
    sol = calc_solar_spectrum(wn,[labpath,'/../',solar_line_file]);
    sol = sol./max(sol); % normalize (ftir spectrum not in physical units)
    
    % lagrange matrix
    wns = (wn-mean(wn))./std(wn);
    wnref = wns([1,fix(end/2),end]);
    L = lagrange(wnref,wns);
    
    disp('Estimating wavelength shift..')
    % wavelength shift
    wn_shift = calc_wn_shift(geo,wn,gasvec,cros,refe,sol,mindep,L);
    
    disp('Estimating solar shift..')
    % solar wl shift
    sol_shift = calc_sol_shift(mfile,wn_shift,labpath,solar_line_file,sol_shift_wn);
    sol = interp1(wn+sol_shift,sol,wn,'linear','extrap');
    
    % default value for offset 
    offset = 0;
    
    % remove some edges of the fitting window
    ncut = 16;
    
    % input for residual / jacobian calculation
    varargin = create_varargin(wn,gasvec,cros,refe,invgas,sol,wn_shift,noise,L,geo,offset,ncut);
    
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
    %nalt = length(geo.center_alts);
    %for n=1:ninvgas
        %P(n) = {ones(nalt,1)};
    %end
    %theta2 = [log(theta(1:ninvgas)); theta(ninvgas+1:end)];
    %[~,~,K1] = jacfun_dr(theta2,ones(ninvgas,1),P,varargin);
    %[out(jj).scaling_A_alpha,out(jj).scaling_A_layer,out(jj).scaling_A_column] = avek_scale(K1,x0,geo.los_lens,varargin{11},geo.altgrid,ncut);
    
    % save results for output
    out(jj).geo = geo;
    out(jj).scaling_factors = theta;
    out(jj).scaling_residual = r;
    out(jj).wn = wn(ncut:end-ncut);
    out(jj).t = refe(ncut:end-ncut);
    out(jj).sol = sol(ncut:end-ncut);
    
    %% ----------------------------------
    %% prior covariance reduction with LM
    %% ----------------------------------
    
    % number of parameters
    d = [k ones(1,ninvgas-1)];
    
    % truncate prior covariance 
    [P, C, Q] = reduce_dim(invgas,d,geo.center_alts);
    
    % initial values
    theta0 = [zeros(1,sum(d)) theta(ninvgas+1:end)'];
    if (fixoff.scale & not(fixoff.dr))
        theta0 = [theta0 offset];
    elseif (not(fixoff.scale) & fixoff.dr)
        theta0 = theta0(1:end-1);
    end
    
    npar = length(theta0);
    
    % jacobian and residual functions
    jacfuni = @(theta0) jacfun_dr(theta0,d,P,varargin);
    resfuni = @(theta0) resfun_dr(theta0,d,P,varargin);
    
    % LM-fit of alpha parameters 
    [theta2,cmat2,rss2,r2] = levmar(resfuni,jacfuni,theta0);
    
    % update error estimate and retrieve again 
    noise = noise*sqrt(rss2);
    varargin = create_varargin(wn,gasvec,cros,refe,invgas,sol,wn_shift,noise,L,geo,offset,ncut);
    jacfuni = @(theta0) jacfun_dr(theta0,d,P,varargin);
    resfuni = @(theta0) resfun_dr(theta0,d,P,varargin);
    [theta2,cmat2,rss2,r2] = levmar(resfuni,jacfuni,theta0);
    
    % averaging kernel
    [~,~,J] = jacfuni(theta2);
    [out(jj).A_alpha,out(jj).A_layer,out(jj).A_column] = avek_dr(J,P{1},theta2(1:d(1)),x0,geo.los_lens,varargin{11},geo.altgrid,ncut);
    
    % retrieved profiles etc. for output
    out(jj).dr_lm_atmos = redu2full(theta2,d,P,invgas,geo.layer_dens);
    out(jj).dr_lm_theta = theta2;
    out(jj).dr_lm_residual = r2;
    out(jj).dr_lm_cmat = cmat2;
    out(jj).dr_pri_C = C;
    out(jj).dr_lm_P = P;
    out(jj).dr_k = k;

    jj = jj + 1;
    
end