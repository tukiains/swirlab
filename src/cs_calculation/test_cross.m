
% Example script for creating Voigt line shapes
% Note: this is terribly slow at the moment and is
% only practical with narrow wavenumber range

clear all
close all

% ch4 window
wnrange = [5996 6008];

% trace gases in this window
gasvec = {'ch4','h2o'};

% what days we create
mdates = {'20131022','20140319'};

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