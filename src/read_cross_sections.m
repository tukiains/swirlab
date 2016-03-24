function [wn,cross,alt]=read_cross_sections(gasvec,wnrange,labpath,voigtpath)

info = load([voigtpath,'voigt_info']);

ind = [];
for n=1:size(wnrange,1)
    ind = cat(2,ind,find(info.wn>=min(wnrange(n,1:2)) & info.wn<=max(wnrange(n,1:2))));
end

wn = info.wn(ind);
wn = wn(:);
alt = info.alt;
clear info

disp('Reading voigt line shapes...')

for n=1:length(gasvec)
    gas = char(gasvec(n));
    disp(gas)
    load([voigtpath,'voigt_shape_',gas]);
    cross(n,:,:) = CS(:,ind);
    clear CS
end

disp('..done')

