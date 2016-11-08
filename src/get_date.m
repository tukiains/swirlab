function d = get_date(fname)
% date = get_date(file_name)

ind = strfind(fname,'so');
ind = ind(end);

d = fname(ind+2:ind+9);