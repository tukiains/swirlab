function vmat = create_varargin(wn,gasvec,cros,refe,invgas,sol,wn_shift,noise,L,geo,offset,ncut,lis)

vmat{1} = wn;
vmat{2} = gasvec;
vmat{3} = cros;
vmat{4} = refe;
vmat{5} = invgas;
vmat{6} = sol;
vmat{7} = wn_shift;
vmat{8} = noise;
vmat{9} = L;
vmat{10} = geo;
vmat{11} = weight_term(sol,noise);
vmat{12} = offset;
vmat{13} = ncut;
vmat{14} = lis;