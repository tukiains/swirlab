function data = read_hitran_cross_section(wnrange, fname)
% data = read_hitran_cross_section(wnrange, filename)

wls = [min(wnrange) max(wnrange)];

load(fname);

wn = data.WaveNumber;

ind = find(wn>min(wnrange) & wn<max(wnrange));

fnames = fieldnames(data);
if (isempty(ind)~=1)
    for n=1:length(fnames)
        temp = data.(char(fnames(n)));
        data.(char(fnames(n))) = temp(ind);
    end
else
    for n=1:length(fnames)        
        data.(char(fnames(n))) = [];
    end
end



