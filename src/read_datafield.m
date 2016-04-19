function data = read_datafield(fname,nskip)
% "generic" reader for any text file that has 
% arbitrary header and data..

fid = fopen(fname);
s = textscan(fid, '%s', 'delimiter', '\n');
fclose(fid);
s = s{1};  % must be cell

% this is slow:
%for n=1:length(s)
%    c = split_string(cell2mat(s(n)));
%    a(n) = sum(cellfun(@(x) length(str2num(x)),c)); % number of numerical values in row
%end

% number of numerical values in each line
% this only works for lines that have only numbers
% should code better version in future
a = cellfun(@(x) length(str2num(x)),s);

% what we have most are data
ind = find(a==mode(a(a>0)));

% read data 
data = cell2mat(cellfun(@str2num,s(ind(nskip:end)),'UniformOutput',false));


