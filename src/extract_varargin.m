function [wn,gasvec,cros,refe,invgas,sol,wn_shift,noise,L,geo] = extract_varargin(a)

a = a{:};

wn = cell2mat(a(1));
gasvec = a(2); gasvec = gasvec{1};
cros = cell2mat(a(3));
refe = cell2mat(a(4));
invgas = a(5); invgas = invgas{1};
sol = cell2mat(a(6));
wn_shift = cell2mat(a(7));
noise = cell2mat(a(8));
L = cell2mat(a(9));
geo = a(10); geo = geo{1};


