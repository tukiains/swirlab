function los_dens_out = redu2full(params,d,P,invgas,los_dens,air)
% los_dens_out = redu2full(params,d,P,invgas,los_dens)

% initialize
los_dens_out = los_dens;

% change retrieved gases
for n = 1:length(invgas)    
    
    reduparams = params(sum(d(1:(n-1)))+1:sum(d(1:n))); 
    
    % Reduced profile for a gas
    dprof = P{n}*reduparams;
    
    gas_old = los_dens.(char(invgas(n))); % in molec
    gas_ppb = gas_old./air*1e9;           % in ppb
    gas_ppb_new = gas_ppb + dprof';       % change 
    gas_new = gas_ppb_new.*air/1e9;       % back to molec

    los_dens_out.(char(invgas(n))) = gas_new;

end