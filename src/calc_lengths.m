function lens = calc_lengths(point,vec,re,altgrid)
% lens = calc_lengths(point,vec,re,altgrid)
%
% algtrid: egdes of the layers (km)

nlay = length(altgrid)-1;

xyz = zeros(nlay,3);
lens = zeros(nlay,1);

for n = 1:nlay
    [~, xyz(n,:)] = crossing_point(point,vec,re,altgrid(n+1));
end
lens(1) = norm(point-xyz(1,:));
for n = 2:nlay
    lens(n) = norm(xyz(n-1,:)-xyz(n,:));
end

