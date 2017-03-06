function [wn,gasvec,cros,refe,invgas,sol,wn_shift,noise,L,geo,err,offset,ncut,lis] = extract_varargin(a)

a = a{:};

wn = a{1};
gasvec = a{2};
cros = a{3};
refe = a{4};
invgas = a{5};
sol = a{6};
wn_shift = a{7};
noise = a{8};
L = a{9};
geo = a{10};
err = a{11};
offset = a{12};
ncut = a{13};
lis = a{14};

