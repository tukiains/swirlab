function [P, C, Q] = reduce_dim(invgas,d,alt)
% [P, C, Q] = reduce_dim(invgas,d,alt)

% Dimensionality reduction for the prior covariance matrix

for n=1:length(invgas)
    
    if (strcmp(invgas(n),'ch4')==1) 
        l = 12; % correlation length
        mu1 = 28; 
        mu2 = 5;  
        sigma1 = 11; 
        sigma2 = 6;
        d1 = 0.4;
        d2 = 0.01;        
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

