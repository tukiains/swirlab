function [P, C, Q] = reduce_dim(invgas,d,alt,ace)
% [P, C, Q] = reduce_dim(invgas,d,alt)

% Dimensionality reduction for the prior covariance matrix
% For LIS use ACE data

for n=1:length(invgas)

    % new covariance marix for ch4 in ppb:
    
    if (strcmp(invgas(n),'ch4')==1) 
        if ace
            % ace profiles over sodankyla
            labpath = fileparts(which('calc_direct_geo.m'));
            ace = load([labpath, '/../input_data/','ace_prior']);
            ace_profs = interp1(ace.alt,ace.acei,alt,'linear','extrap');

            % create prior covariance
            Cor = double(cov(ace_profs'));
        else
            l = 12;
            mu1 = 25;
            mu2 = 5;
            sigma1 = 10;
            sigma2 = 5;
            d1 = 300;
            d2 = 30;
            Cor = create_cov(alt,mu1,mu2,sigma1,sigma2,d1,d2,l);
        end
    else
        Cor = diag(ones(length(alt),1))/1000;
    end

    C{n} = Cor;

    [u,s,v] = svd(Cor);
    Q{n} = diag(s);
    SS = sqrt(Q{n});
    P{n} = u(:,1:d(n))*diag(SS(1:d(n)));
end

