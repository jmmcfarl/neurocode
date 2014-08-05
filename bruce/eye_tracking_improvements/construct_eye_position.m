function [fin_tot_corr,fin_tot_std] = construct_eye_position(fix_post_mean,fix_post_std,drift_post_mean,...
    drift_post_std,fix_ids,trial_start_inds,trial_end_inds,sac_shift)

NT = length(fix_ids);

fin_fix_corr = nan(NT,1);
fin_fix_std = nan(NT,1);
fin_fix_corr(~isnan(fix_ids)) = fix_post_mean(fix_ids(~isnan(fix_ids)));
fin_fix_corr = interp1(find(~isnan(fix_ids)),fin_fix_corr(~isnan(fix_ids)),1:NT);
fin_fix_std(~isnan(fix_ids)) = fix_post_std(fix_ids(~isnan(fix_ids)));
fin_fix_std = interp1(find(~isnan(fix_ids)),fin_fix_std(~isnan(fix_ids)),1:NT);

fin_drift_corr = drift_post_mean;
fin_drift_std = drift_post_std;
fin_drift_corr(isnan(fin_drift_corr)) = 0;
fin_drift_std(isnan(fin_drift_std)) = 0;

for ii = 1:length(trial_start_inds)
    cur_inds = trial_start_inds(ii):trial_end_inds(ii);
    fin_drift_corr(cur_inds(1:end-sac_shift)) = fin_drift_corr(cur_inds(sac_shift+1:end));
    fin_drift_std(cur_inds(1:end-sac_shift)) = fin_drift_std(cur_inds(sac_shift+1:end));
end
fin_drift_corr = interp1(find(~isnan(fix_ids)),fin_drift_corr(~isnan(fix_ids)),1:NT);
fin_drift_std = interp1(find(~isnan(fix_ids)),fin_drift_std(~isnan(fix_ids)),1:NT);

fin_tot_corr = fin_fix_corr + fin_drift_corr;
fin_tot_std = sqrt(fin_fix_std.^2 + fin_drift_std.^2);

%in case you start or within a saccade
fin_tot_corr(isnan(fin_tot_corr)) = 0;
fin_tot_std(isnan(fin_tot_std)) = nan;
