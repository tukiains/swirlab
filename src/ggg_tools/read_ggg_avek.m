function [name,ak,level] = read_ggg_avek(fpath)
% [name,ak,level] = read_ggg_avek(fpath)

mfiles = dir([fpath,'/*.aks']);

for n=1:length(mfiles)
    name(n,:) = mfiles(n).name(1:end-4);
    mfile(n,:) = [fpath mfiles(n).name];

    fid = fopen(mfile(n,:));
    data = textscan(fid,'%d32%f64%f64','headerlines',3);
    fclose(fid);

    ak(n,:) = data{2};    

end

level = data{1};





