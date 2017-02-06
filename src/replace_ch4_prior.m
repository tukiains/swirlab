function geo = replace_ch4_prior(afile,geo)

date = get_date(afile);
ind = strfind(afile,'input_data');

month = str2num(date(5:6));
if (month>=5 & month<=8)
    season = '_summer';
elseif (month>=1 & month<=3)    
    season = '_vortex';
else
    season = '';
end
ch4_prior = load([afile(1:ind+10),'combined_ace_tccon_prior',season,'.mat']);

dd = datenum(str2num(date(1:4)),str2num(date(5:6)),str2num(date(7:8)));
ind = find(ch4_prior.time == dd);

if (ind>0)
    p = ch4_prior.prior(:,ind);
    p2 = interp1(ch4_prior.alt,p,geo.center_alts,'linear','extrap');
    p2 = p2.*geo.air/1e9;
    geo.layer_dens.ch4 = p2;
else
    disp('ERROR: cant find new prior profile for this day!')
    keyboard
end