function r = resfun_dr(theta,d,P,varargin)
% r = resfun_dr(theta,d,P,varargin)
%
% residual function for dimension reduction

[wn,gasvec,cros,refe,invgas,sol,wn_shift,noise,L,geo,err,offset,ncut] = extract_varargin(varargin);

[p1,p2,p3,offset] = fetch_params(theta,invgas,offset,d);

% new profile
dens = redu2full(theta,d,P,invgas,geo.layer_dens);

% transmission
t = calc_direct_radiance(dens,geo.los_lens,gasvec,cros,sol,p1,p2,p3,offset,L);

% convolution
tc = conv_spectrum(wn,t);

% wn shift
tc2 = interp1(wn,tc,wn+wn_shift,'linear','extrap');

% residual 
r = (tc2(:)-refe(:))./err;

% remove edges
r = r(ncut:end-ncut);

% alpha-parameters should be added to residual 
r = [r; theta(1:sum(d))];
