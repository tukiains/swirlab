function show_lm(alt,x0,air,C,theta,cmat,k)
% show_lm(alt,x0,air,C,theta,cmat,k)

a = mvnorr(10000,theta(1:k),cmat(1:k,1:k));

pa = C*a';

profs = bsxfun(@times,exp(pa),x0(:)); 

mix = bsxfun(@rdivide,profs,air(:));

plot_curtain(alt,plims(mix',[0.025 0.5 0.975]),[0.53 0.81 0.98],[0.53 0.81 0.98]);
