function geo = calc_direct_lengths(sza,re,altgrid)

% sun direction
vsun = obj_direction(sza,0);

% instrument location
vi = [0 0 re];

% length of the layers
geo.los_lens = calc_lengths(vi,vsun,re,altgrid);

% save stuff to geo structure
geo.vsun = vsun;
geo.altgrid = altgrid;
geo.re = re;





