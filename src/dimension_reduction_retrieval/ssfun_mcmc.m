function rss = ssfun_mcmc(params,data)
% r = ssfun_mcmc(params,data)

% residual
r = resfun_dr(params(:),data.d,data.P,data.varargin);

% mcmc already has prior for the alpha parameters
r = r(1:end-sum(data.d));

rss = sum(r.^2);

