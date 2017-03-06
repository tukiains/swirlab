function los_dens_out = redu2full(params,d,P,O,invgas,los_dens,air,lis)
% los_dens_out = redu2full(params,d,P,O,invgas,los_dens,air,lis)

% initialize
los_dens_out = los_dens;

% change retrieved gases
for n = 1:length(invgas)    
    
    reduparams = params(sum(d(1:(n-1)))+1:sum(d(1:n))); 
    
    % Reduced profile for a gas
    dprof = P{n}*reduparams;
    
    gas_old = los_dens.(char(invgas(n))); % in molec
    gas_ppb = gas_old./air*1e9;           % in ppb

    if (lis & n==1)
        gas_ppb_new = (P{n}*O*gas_ppb')' + dprof';
    else
        gas_ppb_new = gas_ppb + dprof';
    end

    gas_new = gas_ppb_new.*air/1e9;       % back to molec

    los_dens_out.(char(invgas(n))) = gas_new;

end