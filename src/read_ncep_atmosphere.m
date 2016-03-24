function[air,T,P,atmos] = read_ncep_atmosphere(fname,gasvec,x2,altcorr)

fid = fopen(fname);
s = textscan(fid, '%s', 'delimiter', '\n','headerlines',5);

s = s{1};

for n=1:length(s)
    a = str2num(cell2mat(s(n)));
    alt(n) = a(1) - altcorr; % Sodankylä altitude
    T(n)   = a(2);
    P(n)   = a(3)*1013.25;
    air(n) = a(4);
    h2o(n) = a(5)*air(n);
    co2(n) = a(11)*air(n);
    n2o(n) = a(22)*air(n);
    co(n)  = a(30)*air(n);
    ch4(n) = a(36)*air(n);
    o2(n)  = a(40)*air(n);
    hf(n)  = a(58)*air(n);
    hdo(n) = a(125)*air(n);
end

for n=1:length(gasvec)
    atmos(n,:) = eval(char(gasvec(n)));
end

% interpolate to altitudes that we want
for n=1:length(gasvec)
    atmos2(n,:) = exp(interp1(alt,log(atmos(n,:)),x2,'linear','extrap')); 
end
air = exp(interp1(alt,log(air),x2,'linear','extrap'));
T = interp1(alt,T,x2,'linear','extrap');
P = exp(interp1(alt,log(P),x2,'linear','extrap'));
atmos = atmos2;



