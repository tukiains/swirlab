function ils = ftir_ils(freq,L,th)
% ils = ftir_ils(freq,L,th)
% freq: frequency vector (nm)
% L = OPD of the instrument
% th = FOV in radians
% ils: FTS instrument line shape

gridi = min(freq):0.0001:max(freq);

vn = wl2wn(gridi);          % wl -> wn
v0 = wl2wn(median(gridi));  % center wl

vd = (vn-v0);      % difference vector

vd(vd==0) = 1e-12; % zero difference causes NaN

y1 = sin(2*pi.*vd*L)./(2*pi.*vd*L); % sinc 

if (th==0)
    tils = y1;
else
    y2 = zeros(1,length(y1));
    half = 0.5*v0*th^2;
    ind = find(vd>=-half/2 & vd<=half/2);
    y2(ind) = 2/(v0*th^2);           % rect 
    tils = conv(y2,y1,'same');
end

ils = interp1(gridi,tils,freq,'linear','extrap');

ils = ils./sum(ils);






