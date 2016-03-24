function los_dens_out = redu2full(params,d,P,invgas,los_dens)
% los_dens_out = redu2full(params,d,P,invgas,los_dens)

for n = 1:length(invgas)    
    reduparams = params(sum(d(1:(n-1)))+1:sum(d(1:n))); 
    % Reduced profile for a gas
    dprof = P{n}*reduparams;
    los_dens_out.(char(invgas(n))) = exp(dprof'+log(los_dens.(char(invgas(n))))); % log-normal case
end