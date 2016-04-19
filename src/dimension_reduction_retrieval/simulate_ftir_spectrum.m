function [refe sol_shift wn_shift simuatmos] = simulate_ftir_spectrum(c_wn,cros_o,c_alt,wn,gasvec,afile,sza,L,sol,noise,varargin)

% a fine grid
alt = create_layering(70,100,1.001);

[geo,cros] = calc_direct_geo(c_wn,cros_o,c_alt,wn,gasvec,afile,sza,alt);

p1 = 0.4;
p2 = 0.5;
p3 = 0.5;
offset = 1e-4;

% aircore ch4 profile?
if (nargin>10)

    [~,~,ac_ch4,~,~,~,~,ac_alt] = read_aircore_sounding(varargin{1});
    ac_ch4 = ac_ch4/1e9; % ppb -> mole fraction
    aci = extrapolate_ac(ac_ch4,ac_alt,geo.air,geo.center_alts);
    geo.layer_dens.ch4 = aci;

    % pertubate?
    if (nargin>11)
        per_alt = varargin{2};
        per = varargin{3};        

        [~,ind] = min(abs(per_alt-geo.center_alts));        

        %dx = aci(ind)*per;
        dx = per;
        aci(ind) = aci(ind) + dx;
        
        simuatmos.dx = dx;
        
    end

    geo.layer_dens.ch4 = aci;
    
    
end

simuatmos.ch4 = aci;
simuatmos.alt = geo.center_alts;

% simulate spectrum
refe = calc_direct_radiance(geo.layer_dens,geo.los_lens,gasvec, ...
                            cros,sol,p1,p2,p3,offset,L);

% convolute it
refe = conv_spectrum(wn,refe);

% add noise
rng(666);
refe = refe + noise*randn(length(wn),1);

sol_shift = 0;
wn_shift = 0;