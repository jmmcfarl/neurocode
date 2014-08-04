function filter_contrib = NMMget_filter_contributions(nim, Xstims, stim_filters)
%
% Usage: filter_contrib = NMMget_filter_contributions(nim, Xstims, <stim_filters>)

if nargin < 3
    stim_filters = 1:length(nim.mods);
end

[LL, penLL, pred_rate, G, gint, fgint] = NMMmodel_eval( nim, [], Xstims);

mod_signs = [nim.mods(:).sign];
filtOutMags = std(fgint(:,stim_filters));
tot_stim_filt_out = sum(bsxfun(@times,fgint(:,stim_filters),mod_signs(stim_filters)),2);
totOutMag = std(tot_stim_filt_out);

filter_contrib = filtOutMags/totOutMag;
