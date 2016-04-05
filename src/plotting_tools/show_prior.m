function out = show_prior(alt,x0,air,C)
% show_prior(alt,x0,air,C)
%

a = randn(size(C,2),10000);
pa = C*a;
profs = bsxfun(@times,exp(pa),x0(:)); 
mix = bsxfun(@rdivide,profs,air(:));

out = plot_curtain(alt,plims(mix',[0.025 0.5 0.975]),[.5 .5 .5]);
