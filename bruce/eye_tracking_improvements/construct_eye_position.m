function [fin_tot_corr,fin_tot_std] = construct_eye_position(fix_post_mean,fix_post_std,drift_post_mean,...
    drift_post_std,fix_ids,trial_start_inds,trial_end_inds,neural_delay)
%[fin_tot_corr,fin_tot_std] = construct_eye_position(fix_post_mean,fix_post_std,drift_post_mean,drift_post_std,fix_ids,trial_start_inds,trial_end_inds,neural_delay)
% build mean estimated eye position and its SD from estimates of within-fixation average and drift-corrections
% INPUTS: 
%     fix_pos_mean: [n_fix x 1] vector of posterior mean eye positions during each fixation
%     fix_post_std: [n_fix x 1] vector of posterior SD eye positions during each fixation
%     drift_post_mean: [T x 1] vector of posterior mean drift-corrections (T is total number of time points)
%     drift_post_std: [T x 1] vector of posterior SD drift-corrections
%     fix_ids: [T x 1] vector indexing the fixation number
%     trial_start_inds: [n_trials x 1] vector of trial start indices
%     trial_end_inds: [n_trials x 1] vector of trial end indices
%     neural_delay: scalar number of time bins that estimated EP signals are delayed relative to actual eye positions due to neural repsonse delays
% OUTPUTS: 
%     fin_tot_corr: [T x 1] vector of estimated eye positions (in real time)
%     fin_tot_std: corresponding vector of uncertainties (posterior SD)


NT = length(fix_ids);

fin_fix_corr = nan(NT,1);
fin_fix_std = nan(NT,1);

%back-project data to non-delayed time, and interpolate through saccades/blinks
fin_fix_corr(~isnan(fix_ids)) = fix_post_mean(fix_ids(~isnan(fix_ids)));
fin_fix_corr = interp1(find(~isnan(fix_ids)),fin_fix_corr(~isnan(fix_ids)),1:NT);
fin_fix_std(~isnan(fix_ids)) = fix_post_std(fix_ids(~isnan(fix_ids)));
fin_fix_std = interp1(find(~isnan(fix_ids)),fin_fix_std(~isnan(fix_ids)),1:NT);

fin_drift_corr = drift_post_mean;
fin_drift_std = drift_post_std;
fin_drift_corr(isnan(fin_drift_corr)) = 0;
fin_drift_std(isnan(fin_drift_std)) = 0;
for ii = 1:length(trial_start_inds) %for each trial, back-project inferred drift, and align to start of trial
    cur_inds = trial_start_inds(ii):trial_end_inds(ii);
    fin_drift_corr(cur_inds(1:end-neural_delay)) = fin_drift_corr(cur_inds(neural_delay+1:end));
    fin_drift_std(cur_inds(1:end-neural_delay)) = fin_drift_std(cur_inds(neural_delay+1:end));
end
fin_drift_corr = interp1(find(~isnan(fix_ids)),fin_drift_corr(~isnan(fix_ids)),1:NT); %interpolate through saccades
fin_drift_std = interp1(find(~isnan(fix_ids)),fin_drift_std(~isnan(fix_ids)),1:NT);

%add fixation and drift corrections
fin_tot_corr = fin_fix_corr + fin_drift_corr;
fin_tot_std = sqrt(fin_fix_std.^2 + fin_drift_std.^2); %variances add

%in case you start or within a saccade
fin_tot_corr(isnan(fin_tot_corr)) = 0;
fin_tot_std(isnan(fin_tot_std)) = nan;
