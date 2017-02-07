function [P, C, Q] = reduce_dim(invgas,d,alt)
% [P, C, Q] = reduce_dim(invgas,d,alt)

% Dimensionality reduction for the prior covariance matrix

for n=1:length(invgas)

    % new covariance marix for ch4 in ppb:
    
    if (strcmp(invgas(n),'ch4')==1) 
        l = 12;
        mu1 = 30;
        mu2 = 0;
        sigma1 = 20;
        sigma2 = 0;
        d1 = 2e2;
        d2 = 0;
        Cor = create_cov(alt,mu1,mu2,sigma1,sigma2,d1,d2,l);
    else
        Cor = diag(ones(length(alt),1))/1000;
    end

    C{n} = Cor;

    [u,s,v] = svd(Cor);
    Q{n} = diag(s);
    SS = sqrt(Q{n});
    P{n} = u(:,1:d(n))*diag(SS(1:d(n)));
 
end

