function r = resfun(theta,varargin)
% r=resfuni(theta,varargin)

[wn,gasvec,cros,refe,invgas,sol,wn_shift,noise,L,geo] = extract_varargin(varargin);

p1 = theta(end-3);
p2 = theta(end-2);
p3 = theta(end-1);
offset = theta(end);

% scale densities:
dens = geo.layer_dens;
for n=1:length(invgas)
    dens.(char(invgas(n))) = dens.(char(invgas(n)))*theta(n);
end

% transmission
[t,~,~] = calc_direct_radiance(dens,geo.los_lens,wn,gasvec,cros,sol,p1,p2,p3,offset,L);

% convolution
tc = conv_spectrum(wn,t);

% wn shift
tc2 = interp1(wn,tc,wn+wn_shift,'linear','extrap');

% error
err = weight_term(sol,noise);

% residual
r = (tc2(:)-refe(:))./err;

% ignore edges
r = r(15:end-15);


