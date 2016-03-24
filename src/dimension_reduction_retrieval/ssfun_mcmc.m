function r = ssfun_mcmc(params,data)
% r = ssfun_mcmc(params,data)

% residual
r = resfun_dr(params(:),data.d,data.P,data.varargin);

r = r';
