
clear all
close all

% set path for absorption coeffs:
voigt_path = '/home/tukiains/Dropbox/voigt_shapes/';

[pathstr,name] = fileparts(which('get_ftir_files.m'));

zenlim = 82;
usesimu = false;
lm_only = false;
lis = false;
k = 3;
% fixed offset?
fixo.scale = false; 
fixo.dr = false;
fixo.mcmc = false;

% all dates
datez = {'20130903','20130905','20131022','20140319','20140409','20140508','20140714','20140715','20140716','20140902','20141105'};

for m = 1:length(datez)
    
    mdate = cell2mat(datez(m))
    
    % all files we have
    %prefix = [pathstr,'/../input_data/ftir_spectra/'];
    prefix = ['/home/tukiains/Documents/ggg-2014/ggg-stable/i2s/opus-i2s/spectra/',mdate,'/'];
    
    mfiles = dir([prefix,'so*','*.0*']);
    
    g = 1;
    for n = 1:1:length(mfiles)

        fname = [prefix,mfiles(n).name];
        temp = ftir_dimred_mcmc(voigt_path,fname,lm_only,lis,k,fixo,usesimu,zenlim);

        if (temp.sza>0 & temp.sza<zenlim & isfield(temp,'geo'))
            
            % mean posterior profile
            mcmc_mean_prof(g,:) = mean(temp.mcmc_profs);
            
            % mean column and error
            [mcmc_colu(g) mcmc_colu_err(g)] = posterior_column(temp.mcmc_profs,temp.geo.air,temp.geo.center_alts,[0 20]);
            
            % clear too big fields
            fields = {'mcmc_s2chain','mcmc_chain','mcmc_res','mcmc_profs','dr_pri_C','dr_lm_P'};           
            temp = rmfield(temp,fields);
            
            % save more stuff            
            out(g) = temp;
            mfile(g,:) = fname;
            
            g = g+1;

        end
    end
    
    if (g>1)

        save(['/home/tukiains/data/dimension_reduction_results/retrieved_offset/ftir_results_',mdate,'.mat'],'out','mfile','mcmc_colu','mcmc_colu_err','mcmc_mean_prof');

        clear out mfile mcmc_colu mcmc_colu_err mcmc_mean_prof

    end

end