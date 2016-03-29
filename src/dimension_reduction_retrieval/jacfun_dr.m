function X2 = jacfun_dr(theta,d,P,varargin)
% X2 = jacfun_dr(theta,P,d,varargin)
%
% Jacobian function for dimension reduction

[wn,gasvec,cros,refe,invgas,sol,wn_shift,noise,L,geo] = extract_varargin(varargin);

p1 = theta(end-3);
p2 = theta(end-2);
p3 = theta(end-1);
offset = theta(end);

% new profile
dens = redu2full(theta,d,P,invgas,geo.layer_dens);

% Jacobian
[~,K,K2] = calc_direct_radiance(dens,geo.los_lens,wn,gasvec,cros,sol,p1,p2,p3,offset,L);

% error term
err = weight_term(sol,noise);

% gases
X2 = [];
for n = 1:length(invgas)
    ind = find(ismember(gasvec,invgas(n))==1);
    X = squeeze(K2(ind,:,:))';
    prof = geo.layer_dens.(char(invgas(n))).*geo.los_lens';
    for m=1:size(X,2)
        X(:,m) = conv_spectrum(wn,X(:,m));
        X(:,m) = interp1(wn,X(:,m),wn+wn_shift,'linear','extrap')./err*prof(m);
    end
    X2 = [X2 X*P{n}];
end

% base line + offset
z=1;
for n = sum(d)+1:sum(d)+4
    X2(:,n) = conv_spectrum(wn,K(:,length(invgas)+z));
    X2(:,n) = interp1(wn,X2(:,n),wn+wn_shift,'linear','extrap')./err;
    z=z+1;
end

% remove edges
X2 = X2(15:end-15,:);

% regularization
I = eye(sum(d));
I = [I zeros(sum(d),4)];
X2 = [X2; I];



