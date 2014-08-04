function [fo,COMs,peak_locs] = get_filter_coms_1d(fo)

k_mat = get_k_mat(fo);
[klen,nfilts] = size(k_mat);
SDIM = fo.stim_params.sdim;
flen = klen/SDIM;
x = (1:SDIM)';
COMs = nan(nfilts,1);
peak_locs = nan(nfilts,1);
for i = 1:nfilts
   cur_filt_mat = reshape(fo.mods(i).k,flen,SDIM);
   temp_var_dist = var(cur_filt_mat);
   temp_var_dist = temp_var_dist/sum(temp_var_dist);
   [~,peak_locs(i)] = max(temp_var_dist);
   COMs(i) = temp_var_dist*x;
   fo.mods(i).filt_com = COMs(i);
end