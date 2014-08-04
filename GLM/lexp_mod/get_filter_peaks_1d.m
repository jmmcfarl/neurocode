function [fo,Peaklocs] = get_filter_peaks_1d(fo)

pix_mat = get_pix_mat(fo);
[klen,nfilts] = size(pix_mat);
SDIM = fo.mods(1).SDIM;
flen = klen/SDIM;
x = (1:SDIM)';
Peaklocs = nan(nfilts,1);
for i = 1:nfilts
    cur_filt_mat = reshape(fo.mods(i).pix,flen,SDIM);
    temp_var_dist = var(cur_filt_mat);
    [~,Peaklocs(i)] = max(temp_var_dist);
    fo.mods(i).filt_com = Peaklocs(i);
end