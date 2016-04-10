function [A_alpha,A_layer,A_column] = avek(K,P,theta,x0,los_len,err,alt)
% K: Jacobian
% P: projection matrix
% x0: prior mean
% los_len: lengths of layers
% err: measurement error vector
% alt: altitude grid

x = exp(P*theta(:)).*x0(:); % current profile

Kv = K*diag(los_len); % Jacobian in vertical direction

S = diag(1./err(16:end-16).^2); % inverse of error covariance

Kk =  K*diag(los_len.*x)*P;

A_alpha = (Kk'*S*Kk+eye(length(theta)))\(Kk'*S*Kv); % AK of alpha

A_layer = diag(x)*P*A_alpha; % layer-wise AK

A_column = A_layer*diff(alt)'; % column-wise AK
    