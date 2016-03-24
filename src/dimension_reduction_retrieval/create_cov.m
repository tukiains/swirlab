function C = create_cov(alt,mu1,mu2,sigma1,sigma2,d1,d2,l)

alt = alt(:);
npar = length(alt);

gfun = @(x,l) exp(-0.5*(x./l).^2);
C2 = zeros(npar,npar);
for i=1:npar
    C2(i,:) = gfun(sqrt((alt-alt(i)).^2),l);
end

% sigma "profile"
sig = d1*gfun(alt-mu1,sigma1) + d2*gfun(alt-mu2,sigma2);
C = C2.*(sig*sig');
