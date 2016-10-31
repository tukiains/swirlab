function calc_voigt_shapes_allisotopes(hitran_folder,wn,gasvec,alt,air,T,P,densmat,output_folder)
% calc_voigt_shapes_allisotopes(wn,gasvec,alt,air,T,P,densmat,output_folder)

P_atm = P/1013.25;

for gasloop=1:length(gasvec)
    
    gas = cell2mat(gasvec(gasloop));
    dens = densmat(gasloop,:);
    
    disp('')
    disp(['Calculating Voigt shapes for: ',gas])
    
    savefile = ['voigt_shape_', gas];

    M = molar_mass(gas);

    hfilez = dir([hitran_folder,'hitran_parameters/',gas,'/*.mat']);
    for n=1:length(hfilez)
        filez(n,:) = [hitran_folder,'hitran_parameters/',gas,'/',hfilez(n).name];
    end
    
    CS(1:length(alt),1:length(wn)) = 0;

    for y=1:size(filez,1)

        fname = filez(y,:);
        ind = strfind(fname,'cs');
        isot = str2num(fname(ind+2:end-4));

        disp(['...Isotopologue: ',int2str(isot)])

        data = read_hitran_cross_section([min(wn) max(wn)],fname);

        c = data.LineIntensity;
        freq = data.WaveNumber;             % 1/cm
        tdep = data.TemperatureDependence;
        s_air = data.AirWidth;
        s_self = data.SelfWidth;
        pres = data.PressureShift;
        El = data.LowerStateEnergy;

        Ps = partial_pressure(P,dens,air,'atm');

        % read Q values of this isotope
        qfilee = [hitran_folder,'q_values/',gas,'/',upper(gas),'_',num2str(isot),'.dat'];

        if (isempty(qfilee)==0)
            
            Q = importdata(qfilee);
            Qt = Q(:,1);
            Qv = Q(:,2);
            Q_ref = interp1(Qt,Qv,296); % reference value at 296 K
            Q = interp1(Qt,Qv,T);       % interpolate to this atmosphere temperature
            
            for m=1:length(alt)
                
                % correct line intensity
                c2 = calc_line_intensity_change(freq,c,El,Q_ref,Q(m),T(m));
            
                % correct the line position 
                freq2 = pressure_shift(freq,pres,P_atm(m));
            
                % line broadening factors
                bD = calc_doppler_broadening(T(m),freq2,M);
                bP = calc_pressure_broadening(T(m),P_atm(m),Ps(m),tdep,s_air,s_self);
                
                cross(1:length(wn)) = 0;                        

                for g = 1:length(c2) % all lines 
                
                    ind = findnearest(freq2(g),wn);
                
                    % all wavenumbers 
                    gridi_ind = 1:length(wn);
                    gridi = wn(gridi_ind);
                    
                    % Voigt line shape
                    y = voigt_shape(freq2(g),c2(g),gridi,bP(g),bD(g));                
                    
                    % sum over lines
                    cross(gridi_ind) = cross(gridi_ind) + y;

                end     
       
                % sum over isotopologues
                CS(m,:) = CS(m,:) + cross;            

            end
            
        else
            disp('missing Q values for this isotope!')
        end        
        
    end

    % save results
    cd(output_folder)
    save(savefile,'CS')
    
end

save('voigt_info','alt','wn')

