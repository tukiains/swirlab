function files = list_files(dir_to_path,filetypes)

files = dir([dir_to_path,filetypes]);
files = struct2cell(files);
files = files(1,:);

% join path and file name
for n = 1:length(files)
    files(n) = {[dir_to_path,files{n}]};
end

