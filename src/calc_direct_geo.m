function [geo, cros] = calc_direct_geo(c_wn,cros,c_alt,wn,gasvec,afile,sza,altgrid)

% earth radius 
re = 6377.64;

sza = deg2rad(sza);

% solve geometry
geo = calc_direct_lengths(sza,re,altgrid);
geo.center_alts = calc_center_alts(altgrid);

% interpolate cros sections
cros = interpolate_cross(c_alt, c_wn, cros, geo.center_alts, wn);

% atmosphere
altcorr = 0.184; % Sodankylä altitude
[air,T,P,atmos] = read_ncep_atmosphere(afile,gasvec,altgrid,altcorr);
for n=1:length(gasvec)    
    geo.layer_dens.(char(gasvec(n))) = (atmos(n,1:end-1)+atmos(n,2:end))/2;
end

% air density
geo.air = exp(interp1(altgrid,log(air),geo.center_alts));

