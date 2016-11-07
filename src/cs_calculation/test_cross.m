
% Example script for creating Voigt line shapes
% Note: this is terribly slow at the moment and
% only practical with a narrow wavenumber range

clear all
close all

% ch4 window
%wnrange = [5996 6008];
%gasvec = {'ch4','h2o'};

% o2 window
wnrange = [7870 7890];
gasvec = {'o2','h2o','co2'};

% what days we create
mdates = {'20140319'};

% set me
output_root_folder = '/home/tukiains/data/voigt_test/';

% set me (so*.mav files)
apriori_file_folder = '/home/tukiains/Documents/MATLAB/swirlab_git/swirlab/input_data/ggg_apriori_files/';

% altitude of the TCCON site
altcorr = 0.184; % Sodankylä

% root path of hitran parameters and q-values
% precomputed for hitran 2012 are: https://www.dropbox.com/sh/nuyj46tdp8skpwu/AAApQ3PjvrhnmC6x13GnIFEPa?dl=0
hitran_folder = '/home/tukiains/data/hitran/matfiles/hitran12/';

% create cross sections
create_cross_sections(hitran_folder,apriori_file_folder,output_root_folder,mdates,gasvec,wnrange,altcorr);