function [X2 X3 K1] = jacfun_dr(theta,d,P,O,varargin)
% [J J2 J_first] = jacfun_dr(theta,P,d,varargin)
%
% J: [nwl k] (reduced)
% J2: [nwl (nalt x ngas)] (full)
% J_first: [nwl nalt] (only the first gas)
%
% Jacobian function for dimension reduction

[wn,gasvec,cros,refe,invgas,sol,wn_shift,noise,L,geo,err,offset,ncut,lis] = extract_varargin(varargin);

[p1,p2,p3,offset] = fetch_params(theta,invgas,offset,d);

% new profile
dens = redu2full(theta,d,P,O,invgas,geo.layer_dens,geo.air,lis);

% Jacobian
[~,K,K2] = calc_direct_radiance(dens,geo.los_lens,gasvec,cros,sol,p1,p2,p3,offset,L);

% gases
X2 = [];
X3 = [];
for n = 1:length(invgas)
    ind = find(ismember(gasvec,invgas(n))==1);
    X = squeeze(K2(ind,:,:))';
    %prof = dens.(char(invgas(n))).*geo.los_lens';
    prof = geo.air .* geo.los_lens' / 1e9; % gaussian ppb covariance
    for m = 1:size(X,2)
        X(:,m) = conv_spectrum(wn,X(:,m));
        X(:,m) = interp1(wn,X(:,m),wn+wn_shift,'linear','extrap')./err*prof(m);
    end
    X2 = [X2 X*P{n}];
    X3 = [X3 X];
end

% base line + offset
z=1;
for n = sum(d)+1:length(theta)
    X2(:,n) = conv_spectrum(wn,K(:,length(gasvec)+z));
    X2(:,n) = interp1(wn,X2(:,n),wn+wn_shift,'linear','extrap')./err;
    z = z+1;
end

% remove edges
X2 = X2(ncut:end-ncut,:);
X3 = X3(ncut:end-ncut,:);

% regularization
I = eye(sum(d));
nextra = length(theta)-sum(d);
I = [I zeros(sum(d),nextra)];
X2 = [X2; I];

% only the first gas
if (nargout==3)
    K1 = squeeze(K2(1,:,:));
    for n = 1:size(K1,1)
        K1(n,:) = conv_spectrum(wn,K1(n,:));
        K1(n,:) = interp1(wn,K1(n,:),wn+wn_shift,'linear','extrap');
    end
    K1 = K1(:,ncut:end-ncut)';
end



