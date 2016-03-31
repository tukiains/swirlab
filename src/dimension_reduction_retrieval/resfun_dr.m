function r = resfun_dr(theta,d,P,varargin)
% r = resfun_dr(theta,d,P,varargin)
%
% residual function for dimension reduction

[wn,gasvec,cros,refe,invgas,sol,wn_shift,noise,L,geo,err] = extract_varargin(varargin);

p1 = theta(end-3);
p2 = theta(end-2);
p3 = theta(end-1);
offset = theta(end);

% new profile
dens = redu2full(theta,d,P,invgas,geo.layer_dens);

% transmission
t = calc_direct_radiance(dens,geo.los_lens,wn,gasvec,cros,sol,p1,p2,p3,offset,L);

% convolution
tc = conv_spectrum(wn,t);

% wn shift
tc2 = interp1(wn,tc,wn+wn_shift,'linear','extrap');

% residual 
r = (tc2(:)-refe(:))./err;

% remove edges
r = r(16:end-16);

% alpha-parameters should be added to residual 
r = [r; theta(1:end-4)];
