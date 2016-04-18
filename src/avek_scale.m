function [A_alpha,A_layer,A_column] = avek(K,x0,los_len,err,alt,ncut)

% K: Jacobian [nwl nalt]
% x0: prior mean
% los_len: lengths of layers
% err: measurement error vector
% alt: altitude grid

% dt/dx
Kv = K*diag(los_len);

Ci = diag(1./err(ncut:end-ncut).^2); % inverse of error covariance

Kt = Kv*x0(:);

A_alpha = (Kt'*Ci*Kt)\(Kt'*Ci*Kv); % AK of theta

A_layer = x0(:)*A_alpha; % layer-wise AK

A_column = diff(alt)*A_layer; % column-wise AK 


