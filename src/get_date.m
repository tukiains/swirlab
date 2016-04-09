function d = get_date(fname)
% date = get_date(file_name)

ind = strfind(fname,'/so');

d = fname(ind+3:ind+10);