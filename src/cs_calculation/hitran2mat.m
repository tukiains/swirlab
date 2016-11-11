
clear all
close all

% list of molecules to process
molecules = {'ch4','h2o'};

% path to HITRAN parameter files *.par
% they must be named co2.par, ch4.par, etc.
fpath = '/home/tukiains/data/hitran/hitran_files/2016/';

% output path
opath = '/home/tukiains/data/hitran/matfiles/hitran16/hitran_parameters/';

% format of the *par file
%s = [2 1 12 10 10 5 5 10 4 8 15 15 15 15 ... rest depends on the file];
s = [2 1 12 10 10 5 5 10 4 8]; % these first fields should be (!) always the same
c = cumsum(s);

for loop=1:length(molecules)
    
    mole = cell2mat(molecules(loop));

    disp(['Processing molecule: ', mole, '...'])

    % create directory for this molecule
    mkdir(opath,mole);

    fname = [fpath,mole,'.par'];
    
    % The isotope numbers below are taken from the HITRAN web site.
    % Instead, in the *par files the isotopes are just numbered 1,2,3,...
    % So be careful that the numbers below agree to the *par files you have created
    switch mole
      case 'ch4'
        isotopes = [211 311 212 312];
      case 'h2o'
        isotopes = [161 181 171];
      case 'hdo'
        isotopes = [162 182 172];
      case 'co2'
        isotopes = [626 636 628 627 638 637 828 827 727 838];
      case 'co'
        isotopes = [26 36 28 27 38 37];
      case 'o2'
        isotopes = [66 68 67];
    end        

    fid = fopen(fname,'r');
    temp = textscan(fid,'%s','delimiter','\n','Whitespace','');
    fclose(fid);

    temp = char(temp{1});
   
    % information we need 
    isotope     = str2num(temp(:,[c(1)+1:c(2)]));   % isotopologue number 
    center_freq = str2num(temp(:,[c(2)+1:c(3)]));   % transition wavenumber (1/cm)
    intens      = str2num(temp(:,[c(3)+1:c(4)]));   % line intensity at 296 K
    a_coeff     = str2num(temp(:,[c(4)+1:c(5)]));   % einstein coeff
    air_width   = str2num(temp(:,[c(5)+1:c(6)]));   % air-broaneded width
    self_width  = str2num(temp(:,[c(6)+1:c(7)]));   % self-broadened width
    lower_e     = str2num(temp(:,[c(7)+1:c(8)]));   % lower-state energy
    temp_dep    = str2num(temp(:,[c(8)+1:c(9)]));   % T-dependency (of air width)
    press_shift = str2num(temp(:,[c(9)+1:c(10)]));  % pressure shift
    
    for g=1:length(isotopes)               
        
        ind = find(isotope==g);
        
        % current isotope
        center_freq2 = center_freq(ind);
        intens2      = intens(ind);        
        a_coeff2     = a_coeff(ind);
        air_width2   = air_width(ind);
        self_width2  = self_width(ind);
        lower_e2     = lower_e(ind);
        temp_dep2    = temp_dep(ind);
        press_shift2 = press_shift(ind);        
    
        % there are a few duplicate lines (should we take them all or remove duplicates?)
        [center_freq2 ind] = unique(center_freq2,'first');
        intens2      = intens2(ind);
        a_coeff2     = a_coeff2(ind);
        air_width2   = air_width2(ind);
        self_width2  = self_width2(ind);
        lower_e2     = lower_e2(ind);
        temp_dep2    = temp_dep2(ind);
        press_shift2 = press_shift2(ind);
        
        data = struct('WaveNumber',center_freq2,...
                      'LineIntensity',intens2,...
                      'ACoefficient', a_coeff2,...
                      'AirWidth',air_width2,...
                      'SelfWidth',self_width2, ...
                      'TemperatureDependence',temp_dep2,...
                      'PressureShift',press_shift2, ...
                      'LowerStateEnergy',lower_e2);
        
        save([opath, mole, '/', upper(mole), 'cs', int2str(isotopes(g))],'data')
        
    end 

    loop = loop+1;
    
end

