function [all_tkerns,avg_Ekern,avg_Ikern] = get_hilbert_tempkerns(stim_mod)

stim_mod_signs = [stim_mod.mods(:).sign];
nmods = length(stim_mod_signs);
stim_dims = stim_mod.stim_params(1).stim_dims;
flen = stim_dims(1);

mod_filts = reshape([stim_mod.mods(:).filtK],[stim_dims(1) stim_dims(2) nmods]);

all_tkerns = nan(flen,nmods);
for jj = 1:nmods
    all_tkerns(:,jj) = squeeze(mean(abs(hilbert(squeeze(mod_filts(:,:,jj)))),2));
end

%rescale based on relative contributions of each filter
all_tkerns = bsxfun(@rdivide,all_tkerns,sqrt(sum(all_tkerns.^2)));
all_tkerns = bsxfun(@times,all_tkerns,stim_mod.rel_filt_weights);

Emods = find(stim_mod_signs == 1);
Imods = find(stim_mod_signs == -1);
avg_Ekern = mean(all_tkerns(:,Emods),2)/mean(stim_mod.rel_filt_weights(Emods));
avg_Ikern = mean(all_tkerns(:,Imods)/2)/mean(stim_mod.rel_filt_weights(Imods));