function [X2 X3] = jacfun_dr(theta,d,P,varargin)
% [J J2] = jacfun_dr(theta,P,d,varargin)
%
% J: [nwl k] (reduced)
% J2: [nwl nalt] (full)
%
% Jacobian function for dimension reduction

[wn,gasvec,cros,refe,invgas,sol,wn_shift,noise,L,geo,err] = extract_varargin(varargin);

p1 = theta(end-3);
p2 = theta(end-2);
p3 = theta(end-1);
offset = theta(end);

% new profile
dens = redu2full(theta,d,P,invgas,geo.layer_dens);

% Jacobian
[~,K,K2] = calc_direct_radiance(dens,geo.los_lens,wn,gasvec,cros,sol,p1,p2,p3,offset,L);

% gases
X2 = [];
X3 = [];
for n = 1:length(invgas)
    ind = find(ismember(gasvec,invgas(n))==1);
    X = squeeze(K2(ind,:,:))';
    prof = dens.(char(invgas(n))).*geo.los_lens';
    for m = 1:size(X,2)
        X(:,m) = conv_spectrum(wn,X(:,m));
        X(:,m) = interp1(wn,X(:,m),wn+wn_shift,'linear','extrap')./err*prof(m);
    end
    X2 = [X2 X*P{n}];
    X3 = [X3 X];
end

% base line + offset
z=1;
for n = sum(d)+1:sum(d)+4
    X2(:,n) = conv_spectrum(wn,K(:,length(gasvec)+z));
    X2(:,n) = interp1(wn,X2(:,n),wn+wn_shift,'linear','extrap')./err;
    z = z+1;
end

% remove edges
X2 = X2(16:end-16,:);
X3 = X3(16:end-16,:);

% regularization
I = eye(sum(d));
I = [I zeros(sum(d),4)];
X2 = [X2; I];





