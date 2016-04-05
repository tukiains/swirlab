function out = plot_curtain(alt,lims,varargin)

if (nargin>2)
    colos = varargin{1};
    if (nargin==4)
        colom = varargin{2};
    end        
else
    colos = [0.6 0.6 0.6];
end

np = size(lims,1);
nn = (np+1)/2;    % median

h = fill_curtain(alt,lims(1,:),lims(2*nn-1,:),colos);

hold on
for k=2:(nn-1)
    h(n) = fill_curtain(alt,lims(k,:),lims(2*nn-k,:),dimc.*0.9.^(k-1));
end

% plot median?
if (nargin==4)
    plot(lims(nn,:),alt,'linestyle','-','color',colom,'linewidth',2,'handlevisibility','off');
end

if (nargout>0)
    out = h;
end