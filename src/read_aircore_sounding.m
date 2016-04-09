function [co2,co2e,ch4,ch4e,co,coe,pres,alt,temp,air] = read_aircore_sounding(fname)
% [co2,co2e,ch4,ch4e,co,coe,pres,alt,temp,air] = read_aircore_sounding(fname)
%

fid = fopen(fname);
s = textscan(fid, '%s', 'delimiter', '\n');
fclose(fid);
s = s{1};

% number of numerical values in each line
a = cellfun(@(x) length(str2num(x)),s);

% what we have most are data
ind = find(a==mode(a(a>0)));

% read data (skip 10 first lines to make sure we dont include any header data)
data = cell2mat(cellfun(@str2num,s(ind(10:end)),'UniformOutput',false));

data(data==99999) = nan;

% data fields
pres = data(:,2);
co2  = data(:,3);
co2e = 2*data(:,6);
ch4  = data(:,4);
ch4e = 2*data(:,7);
co   = data(:,5);
coe  = 2*data(:,8);
alt  = data(:,11)/1000;
temp = data(:,9);

% air from ideal gas law
R = 287.058;
air = pres*10./(R.*temp)*1e6;
