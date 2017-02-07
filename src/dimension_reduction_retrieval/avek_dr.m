function [A_alpha,A_layer,A_column] = avek_dr(K,P,theta,x0,los_len,err,alt,ncut,air)

% K: Jacobian
% P: projection matrix
% x0: prior mean
% los_len: lengths of layers
% err: measurement error vector
% alt: altitude grid
% air: air profile to covert

x = P*theta(:) + x0(:); % current profile

% dt/dx
Kv = K * diag(los_len.*air'/1e9); % -> into ppb space?

S = diag(1./err(ncut:end-ncut).^2); % inverse of error covariance

% dt/da
Ka =  K*diag(los_len.*air'/1e9)*P; % ppb space

A_alpha = (Ka'*S*Ka+eye(length(theta)))\(Ka'*S*Kv); % AK of alpha

A_layer = P*A_alpha; % layer-wise AK

A_column = diff(alt)*A_layer;    

