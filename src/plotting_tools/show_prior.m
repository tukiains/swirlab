function [out mix] = show_prior(alt,x0,air,C)
% show_prior(alt,x0,air,C)
%

a = randn(size(C,2),10000);
pa = C*a;

mix = pa + x0(:)./air'*1e9;

out = plot_curtain(alt,plims(mix',[0.025 0.5 0.975]),[.5 .5 .5]);
