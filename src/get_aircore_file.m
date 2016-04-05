function fname = get_aircore_file(mdate,fpath)
% fname = get_aircore_file(mdate,fpath)
%
% mdate: 'yyyymmdd'
% fpath: path to aircore files
% fname: full path of the aircore file

files = dir([fpath,'*',mdate,'*.naf']);

if (length(files)>0)
    % take first (there can be several)
    fname = [fpath,files(1).name];
end

