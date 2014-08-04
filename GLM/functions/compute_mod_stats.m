function [spatial_profiles, temporal_profiles, weights, mod_type, space_COM, temp_COM] = ...
    compute_mod_stats(fo)
% USAGE: RF = plot2drf(k,flen)
% Plots the RF and gives a matrix containing it back

klen = length(fo.mods(1).k);
sdim = fo.mods(1).SDIM;
flen = klen/sdim;
nmods = length(fo.mods);

spatial_profiles = nan(nmods,sdim);
temporal_profiles = nan(nmods,flen);
weights = nan(nmods,1);
space_COM = nan(nmods,1);
temp_COM = nan(nmods,1);
mod_type = nan(nmods,1);
space_ax = 1:sdim;
temp_ax = 1:flen;


for imod = 1:nmods
    
   cur_RF = reshape(fo.mods(imod).k,flen,sdim);
   spatial_profiles(imod,:) = var(cur_RF);
   temporal_profiles(imod,:) = var(cur_RF,[],2);
   weights(imod) = fo.mods(imod).w;
   space_COM(imod) = sum(space_ax.*spatial_profiles(imod,:))/sum(spatial_profiles(imod,:));
   temp_COM(imod) = sum(temp_ax.*temporal_profiles(imod,:))/sum(temporal_profiles(imod,:));
    
   if mean(fo.mods(imod).nly) > 0
       mod_type(imod) = 1;
   elseif mean(fo.mods(imod).nly) < 0
       mod_type(imod) = -1;
   end
   
end

