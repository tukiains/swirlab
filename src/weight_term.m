function [err] = weight_term(sol,noise)

err = noise*ones(length(sol),1);

% add extra weigth to strong solar lines?
%ind = find(sol<0.94);
%err(ind) = 0.004./sol(ind);



