function plot_curtain(alt,lims,varargin)

if (nargin>2)
    colom = varargin{1};
    colos = varargin{2};
else
    colom = [0.53 0.81 0.98];
    colos = [0.6 0.6 0.6];
end

np = size(lims,1);
nn = (np+1)/2; % median

fill_curtain(alt,lims(1,:),lims(2*nn-1,:),colos);

hold on
for k=2:(nn-1)
    fill_curtain(alt,lims(k,:),lims(2*nn-k,:),dimc.*0.9.^(k-1));
end

if (nargin>2)
    plot(lims(nn,:),alt,'linestyle','-','color',colom,'linewidth',2,'handlevisibility','off');
end