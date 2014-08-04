function [out_of_range] = detect_bad_fixation(corrected_interp_eyevals,all_trialvec,used_inds,par_thresh,orth_thresh,target_coils)

if nargin < 6
    target_coils = [1 0];
end

if target_coils(1) == 1 && target_coils(2) == 0
    corrected_interp_eyevals = corrected_interp_eyevals(:,[1 2]);
elseif target_coils(1) == 0 && target_coils(2) == 1
    corrected_interp_eyevals = corrected_interp_eyevals(:,[3 4]);
elseif target_coils(1) == 1 && target_coils(2) == 1    
    corrected_interp_eyevals = 0.5*corrected_interp_eyevals(:,[3 4])+0.5*corrected_interp_eyevals(:,[1 2]);
else
    error('Invalid target coils');
end
    
%% DETECT TIMES WHEN TRIAL FAILURES OCCUR AND FORCE-END TRIALS
fprintf('Looking for fixation breaks\n');

back_samps = 1;

% out_sapmles = find(abs(corrected_interp_eyevals(:,2)) >= orth_thresh | ...
%     abs(corrected_interp_eyevals(:,1)) >= par_thresh);

un_trials = unique(all_trialvec(used_inds))';
out_of_range = [];
for ii = un_trials
%     fprintf('Trial %d of %d\n',ii,length(un_trials));
    cur_trial_inds = used_inds(all_trialvec(used_inds) == ii);
    if ~isempty(cur_trial_inds)
        cur_out = find(abs(corrected_interp_eyevals(cur_trial_inds,2)) >= orth_thresh | ...
            abs(corrected_interp_eyevals(cur_trial_inds,1)) >= par_thresh,1,'first');
        if ~isempty(cur_out)
            if cur_out > back_samps
                cur_out = cur_out - back_samps;
            end
            cur_out = cur_out:length(cur_trial_inds);
        end
        out_of_range = cat(1,out_of_range,cur_trial_inds(cur_out(:)));
    end
end


