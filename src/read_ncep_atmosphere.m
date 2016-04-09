function[air,T,P,atmos] = read_ncep_atmosphere(fname,gasvec,alt_in,altcorr)

nhead = 7; % do we have to fix this?

fid = fopen(fname);
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

for n=1:length(s)
    
    a = str2num(cell2mat(s(n)));

    alt(n) = a(fn(fields, 'Height')) - altcorr; 
    T(n)   = a(fn(fields, 'Temp'));
    P(n)   = a(fn(fields, 'Pres')) * 1013.25;
    air(n) = a(fn(fields, 'Density'));
    h2o(n) = a(fn(fields, '1h2o'));
    co2(n) = a(fn(fields, '1co2'));
    n2o(n) = a(fn(fields, '1n2o'));
    co(n)  = a(fn(fields, '1co'));
    ch4(n) = a(fn(fields, '1ch4'));
    o2(n)  = a(fn(fields, '1o2'));
    hf(n)  = a(fn(fields, '1hf'));
    hdo(n) = a(fn(fields, '1hdo'));

end

for n=1:length(gasvec)
    prof = eval(char(gasvec(n))) .* air;    
    atmos(n,:) = exp(interp1(alt,log(prof),alt_in,'linear','extrap'));
end

air = exp(interp1(alt,log(air),alt_in,'linear','extrap'));
T = interp1(alt,T,alt_in,'linear','extrap');
P = exp(interp1(alt,log(P),alt_in,'linear','extrap'));



