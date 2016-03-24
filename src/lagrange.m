function L=lagrange(xref,x)
%LAGRANGE calculate Lagrange polynomials
% L=lagrange(xref,x), xref: data points, x evaluation points
% returns matrix L for which y = L*yref, solves interpolating polynomial at x
    
n = length(xref);
m = length(x);
L = ones(m,n);
ii = 1:n;

for i=1:n
    for j = ii(ii~=i)
        L(:,i) = L(:,i).*(x(:)-xref(j))./(xref(i)-xref(j));
    end
end