function calt = calc_center_alts(altgrid);
% calt = calc_center_alts(altgrid);

calt = altgrid(1:(end-1)) + diff(altgrid)/2;
