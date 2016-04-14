function [p1,p2,p3,offset] = fetch_params(theta,invgas,offset,varargin);
% [p1,p2,p3,offset] = fetch_params(theta,invgas,offset,d);

if (nargin==3)
    npar = length(invgas);
elseif (nargin==4)
    d = varargin{1};
    npar = sum(d);
end
    
% polynom terms
p1 = theta(npar+1);
p2 = theta(npar+2);
p3 = theta(npar+3);

% offset term (may be fixed)
if (length(theta)==(npar+4))
    offset = theta(end);
end

