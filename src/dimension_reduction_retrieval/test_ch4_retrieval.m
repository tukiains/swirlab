
clear all
close all

% path for absorption coeffs:
voigt_path = '/home/tukiains/data/voigt_shapes/';

% measurement file:
mfile = '/home/tukiains/Documents/MATLAB/swirlab_git/swirlab/input_data/ftir_spectra/wg20140319saebaa.067';

% do retrieval
out = ftir_dimred_mcmc(voigt_path,mfile);

% show results

