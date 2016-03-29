function convoluted_spec = conv_spectrum(wn,spec)
% convoluted_spectrum = conv_spectrum(wn,spectrum)

persistent f;

if isempty(f) 
    % Sodankylä FTS: OPD 45 cm, FOV 2.3923 mrad
    f = ftir_ils(wn2wl(wn),45,2.3923e-3);
end

convoluted_spec = conv(spec,f,'same');


    
