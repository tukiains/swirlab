function convoluted_spec = conv_spectrum(wn,spec)
% convoluted_spectrum = conv_spectrum(wn,spectrum)

persistent f;

if isempty(f) 
    % Sodankylä FTS: 
    opd = 45;
    fov = 2.3923e-3;
    f = ggg_ils(2,length(wn),median(wn),median(diff(wn)),opd,fov);
end

convoluted_spec = conv(spec,f,'same');


    
