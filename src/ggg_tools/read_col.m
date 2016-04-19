function [name, sza, VSF_ch4, OVC_ch4, VSF_h2o, OVC_h2o] = read_col(colfile,grlfile)
% [name, sza, VSF_ch4, OVC_ch4, VSF_h2o, OVC_h2o] = read_col(colfile,grlfile)
% read ggg.col file

nhead = 24; % do we have to fix this?

fid = fopen(colfile);
% read names of the gases
for n=1:nhead
    fields = fgets(fid);
end
fields = split_string(fields);

s = textscan(fid, '%s', 'delimiter', '\n');
s = s{1};
fclose(fid);

% to find correct field number
fn = @(fields,name) find(ismember(fields,name)==1);

% this loop is slow and stupid
for n=1:length(s)
    a = split_string(cell2mat(s(n)));
    
    b = a(fn(fields, 'VSF_ch4'));
    VSF_ch4(n) = str2num(cell2mat(b));
    
    b = a(fn(fields, 'OVC_ch4'));
    OVC_ch4(n) = str2num(cell2mat(b));
    
    b = a(fn(fields, 'VSF_h2o'));
    VSF_h2o(n) = str2num(cell2mat(b));

    b = a(fn(fields, 'OVC_h2o'));
    OVC_h2o(n) = str2num(cell2mat(b));

    name(n,:) = a(1);
    
end


% info from grl file
[sza, names] = read_zenith_angles(grlfile,90);
sza = sza(1:2:end);

