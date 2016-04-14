function [A_alpha,A_layer,A_column] = avek_dr(K,P,theta,x0,los_len,err,alt,ncut)

% K: Jacobian
% P: projection matrix
% x0: prior mean
% los_len: lengths of layers
% err: measurement error vector
% alt: altitude grid

x = exp(P*theta(:)).*x0(:); % current profile

% dt/dx
Kv = K * diag(los_len);

S = diag(1./err(ncut:end-ncut).^2); % inverse of error covariance

% dt/da
Ka =  K*diag(los_len.*x)*P; 

A_alpha = (Ka'*S*Ka+eye(length(theta)))\(Ka'*S*Kv); % AK of alpha

A_layer = diag(x)*P*A_alpha; % layer-wise AK

A_column = A_layer*diff(alt)'; % column-wise AK
    