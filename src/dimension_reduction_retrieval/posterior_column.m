function [colum colu_err] = posterior_column(profs,air,alt,altlim)

alt2 = [min(altlim):0.1:max(altlim)];

prof = bsxfun(@times,profs,air);
profi = exp(interp1(alt,log(prof)',alt2,'linear','extrap'));

colu = trapz(alt2,profi);
colum = mean(colu);
colu_err = 2*std(colu);

