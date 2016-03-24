function fname = get_ftir_files(mdate,fpath)
% fname = get_aircore_file(mdate,fpath)
%
% mdate: 'yyyymmdd'
% fpath: path to ftir files
% fname: full path of the ftir files

files = dir([fpath,'*',mdate,'*']);

for n=1:length(files)
    fname(n) = {[fpath,files(n).name]};
end
