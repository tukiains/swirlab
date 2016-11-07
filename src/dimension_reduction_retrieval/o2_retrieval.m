function out = o2_retrieval(voigt_root_path,mfile,zenlim)

% date of the measurement
mdate = get_date(mfile);

% rootpath of swirlab
labpath = fileparts(which('calc_direct_geo.m'));

% select window 
[window,wnrange,gasvec,invgas,ninvgas,sol_shift_wn,solar_line_file,mindep] = window_details('o2',1);

% solar zenith angle
out.sza = get_sza_angle(mfile,[labpath,'/../input_data/ggg_runlog_files/so',mdate,'.grl'])

if (isempty(out.sza) | out.sza>zenlim)
    disp('Error: cant find solar zenith angle for this scan (or too high)')
    return
end

% voig path
voigtpath = [voigt_root_path,'voigt_shapes_',mdate,'_', ...
             int2str(min(window)),'-',int2str(max(window)),'/'];

% cross sections 
[c_wn,cros_o,c_alt] = read_cross_sections(gasvec, wnrange, labpath, voigtpath);

% ggg prior atmosphere file
afile = [labpath,'/../input_data/ggg_apriori_files/so',mdate,'.mav'];

% noise of the spectrum 
noise = 4e-4;

% reference (FTIR measurement)
[refe,wn] = read_ftir_spectrum(mfile,wnrange);

if (isempty(refe))
    disp('Error: No spectrum')
    return
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

% wavelength shift
wn_shift = calc_wn_shift(geo,wn,gasvec,cros,refe,sol,mindep,L);

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
theta0 = [ones(1,ninvgas), 0.5, 0.5, 0.5, offset]; % must also retrieve offset

% LM-fit of scaling factors
[theta,cmat,rss,r] = levmar(@resfun,@jacfun,theta0,varargin);

% save results for output
out.geo = geo;
out.scaling_factors = theta;
out.scaling_residual = r;
out.wn = wn(ncut:end-ncut);
out.t = refe(ncut:end-ncut);
out.sol = sol(ncut:end-ncut);

