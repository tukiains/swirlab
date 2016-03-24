function altgrid = create_layering(altmax,N,b)
% altgrid = create_layering(altmax,N,b)
% layring between 0..altmax with N layers
% and assymetry factor b

x1 = 0:N-1;
a1 = altmax * (1-b)/(1-b^N);
altgrid = [0 cumsum(a1*b.^x1)];
