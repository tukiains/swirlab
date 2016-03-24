function cros_i = interpolate_cross(alt,wn,cros,alt2,wn2);
% cros_i = interpolate_cross(alt,wn,cros,alt2,wn2);

for n=1:size(cros,1)

    a = squeeze(cros(n,:,:));
    b = interp2(alt,wn,a',alt2,wn2,'linear',0);
    cros_i(n,:,:) = b';

end
