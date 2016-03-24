function Ko = jacfun(theta,varargin)
% K=jacfun(theta,varargin)

[wn,gasvec,cros,refe,invgas,sol,wn_shift,noise,L,geo] = extract_varargin(varargin);

p1 = theta(end-3);
p2 = theta(end-2);
p3 = theta(end-1);
offset = theta(end);

% number of retrieved gases
ninvgas = length(invgas);

% scale densities:
dens = geo.layer_dens;
for n=1:ninvgas
    dens.(char(invgas(n))) = dens.(char(invgas(n)))*theta(n);
end

% evaluate Jacobian
[~,K,~] = calc_direct_radiance(dens,geo.los_lens,wn,gasvec,cros,sol,p1,p2,p3,offset,L);

% error
err = weight_term(sol,noise);

% only retrieved gases
for n = 1:ninvgas
    ind = find(ismember(gasvec,invgas(n))==1);
    Ko(:,n) = conv_spectrum(wn,K(:,ind));
    Ko(:,n) = interp1(wn,Ko(:,n),wn+wn_shift,'linear','extrap');
    Ko(:,n) = Ko(:,n)./err;
end

% assumes 4 extra retrieved parameters:
for n = ninvgas+1:ninvgas+4
    ind = length(gasvec)+n-ninvgas;
    Ko(:,n) = conv_spectrum(wn,K(:,ind));
    Ko(:,n) = interp1(wn,Ko(:,n),wn + wn_shift, 'linear', 'extrap');
    Ko(:,n) = Ko(:,n)./err;
end

% ignore edges
Ko = Ko(15:end-15,:);

