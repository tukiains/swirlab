function convoluted_spec = conv_spectrum(wn,spec)
% convoluted_spectrum = conv_spectrum(wn,spectrum)

persistent f;

if isempty(f) 
    f = ggg_ils(2,length(wn),median(wn),median(diff(wn)),45,2.3923e-3);
    f = f./sum(f);
end
convoluted_spec = conv(spec,f,'same');


    
