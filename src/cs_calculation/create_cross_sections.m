function create_cross_sections(hitran_folder,apriori_file_folder,output_root_folder,mdates,gasvec,wnrange,altcorr)

labpath = fileparts(which('calc_direct_geo.m'));

% step in Voigt calculation [1/cm]
step = 0.0005;

% wn vector
wn = [min(wnrange)-3:step:max(wnrange)+3];

% layering of the atmosphere 
alt = create_layering(70,70,1.02);

% loop over dates
for n=1:length(mdates)

    mdate = cell2mat(mdates(n))

    afile = [apriori_file_folder,'so',mdate,'.mav'];

    % output folder
    output_folder = [output_root_folder,'voigt_shapes_',mdate,'_',int2str(min(wnrange)),'-',int2str(max(wnrange)),'/'];
    
    % create it
    if (exist(output_folder)~=7)
        system(['mkdir ',output_folder]);
    end

    % read atmosphere
    [air,T,P,atmos] = read_ncep_atmosphere(afile,gasvec,alt,altcorr);

    % calculate cross sections
    calc_voigt_shapes_allisotopes(hitran_folder,wn,gasvec,alt,air,T,P,atmos,output_folder);

end
