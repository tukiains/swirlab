function [t,K,K2] = calc_direct_radiance(layer_dens,los_lens,wn,gasvec,cros,sol,p1,p2,p3,offset,L)

% km to cm
c = 1e5;

% parametrized polynomial and solar sepctrum
base = L*[p1 p2 p3]' .* sol;

% size of the problem
[ngas,nlay,nwl] = size(cros);

% densities of the layers
dens_mat = cell2mat(struct2cell(layer_dens));

% .. as vector
dens_lay = reshape(dens_mat,ngas*nlay,1);

% cross section matrix
cmat = c*reshape(cros,ngas*nlay,nwl);

% density * length
dens_los = reshape(repmat(los_lens,1,ngas)',1,nlay*ngas)'.*dens_lay;

% transmission
t_beam = exp(-cmat'*dens_los);

% add baseline and offset 
t = t_beam .* base + offset;

% derivatives part I: 1 values / wl / gas
for n=1:ngas
    K(:,n) = -c * squeeze(cros(n,:,:))'*(dens_mat(n,:)'.*los_lens) .*t_beam .* base; % gases
end
K = [K diag(t_beam)*diag(sol)*L ones(length(wn),1)]; % base line

% derivatives part II: 1 value / wl / gas / layer
G = reshape(diag(dens_los)*cmat,ngas,nlay,nwl);
for n=1:ngas
    K2(n,:,:) = bsxfun(@times,-squeeze(G(n,:,:)),(t_beam .* base)'); % just gases
end





