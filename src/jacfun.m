function Ko = jacfun(theta,varargin)
% K=jacfun(theta,varargin)

[wn,gasvec,cros,refe,invgas,sol,wn_shift,noise,L,geo,err,offset,ncut] = extract_varargin(varargin);

ninvgas = length(invgas);

[p1,p2,p3,offset] = fetch_params(theta,invgas,offset);

% scale densities:
dens = geo.layer_dens;
for n = 1:ninvgas
    dens.(char(invgas(n))) = dens.(char(invgas(n)))*theta(n);
end

% evaluate Jacobian
[~,K] = calc_direct_radiance(dens,geo.los_lens,gasvec,cros,sol,p1,p2,p3,offset,L);

% retrieved gases
for n = 1:ninvgas
    ind = find(ismember(gasvec,invgas(n))==1);
    Ko(:,n) = conv_spectrum(wn,K(:,ind));
    Ko(:,n) = interp1(wn,Ko(:,n),wn+wn_shift,'linear','extrap');
    Ko(:,n) = Ko(:,n)./err;
end

% polynom temrms and possible the offset term
for n = ninvgas+1:length(theta)
    ind = length(gasvec)+n-ninvgas;
    Ko(:,n) = conv_spectrum(wn,K(:,ind));
    Ko(:,n) = interp1(wn,Ko(:,n),wn + wn_shift, 'linear', 'extrap');
    Ko(:,n) = Ko(:,n)./err;
end

% ignore edges
Ko = Ko(ncut:end-ncut,:);

