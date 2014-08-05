
info_backlag = round(0.3/dt);
info_forlag = round(0.3/dt);
info_slags = -info_backlag:info_forlag;

cur_sac_start_inds = saccade_start_inds(big_sacs);
cur_sac_stop_inds = saccade_stop_inds(big_sacs);

badsacs = find(cur_sac_start_inds(2:end) - cur_sac_stop_inds(1:end-1) < 50);
cur_sac_start_inds(badsacs+1) = [];
cur_sac_stop_inds(badsacs+1) = [];

rand_jit = randi(100,length(cur_sac_start_inds),1);
rcur_sac_start_inds = cur_sac_start_inds + rand_jit;
rcur_sac_stop_inds = cur_sac_stop_inds + rand_jit;
% rcur_sac_stop_inds = rcur_sac_start_inds;
bad = find(rcur_sac_stop_inds > NT);
rcur_sac_start_inds(bad) = [];
rcur_sac_stop_inds(bad) = [];

saccade_stop_trial_inds = all_trialvec(used_inds(cur_sac_stop_inds));
rsaccade_stop_trial_inds = all_trialvec(used_inds(rcur_sac_stop_inds));

Xsac_end = zeros(NT,length(info_slags));
Xsac_rend = zeros(NT,length(info_slags));
for ii = 1:length(info_slags)
    cur_sac_target = cur_sac_stop_inds + info_slags(ii);
    uu = find(cur_sac_target > 1 & cur_sac_target < NT);
    cur_sac_target = cur_sac_target(uu);
    cur_sac_target(all_trialvec(used_inds(cur_sac_target)) ~= saccade_stop_trial_inds(uu)) = [];
    Xsac_end(cur_sac_target,ii) = 1;
    
    cur_sac_target = rcur_sac_stop_inds + info_slags(ii);
    uu = find(cur_sac_target > 1 & cur_sac_target < NT);
    cur_sac_target = cur_sac_target(uu);
    cur_sac_target(all_trialvec(used_inds(cur_sac_target)) ~= rsaccade_stop_trial_inds(uu)) = [];
    Xsac_rend(cur_sac_target,ii) = 1;
end
Xsac_end = Xsac_end(cc_uinds,:);
Xsac_rend = Xsac_rend(cc_uinds,:);

Xsac_start = zeros(NT,length(info_slags));
Xsac_rstart = zeros(NT,length(info_slags));
for ii = 1:length(info_slags)
    cur_sac_target = cur_sac_start_inds + info_slags(ii);
    uu = find(cur_sac_target > 1 & cur_sac_target < NT);
    cur_sac_target = cur_sac_target(uu);
    cur_sac_target(all_trialvec(used_inds(cur_sac_target)) ~= saccade_stop_trial_inds(uu)) = [];
    Xsac_start(cur_sac_target,ii) = 1;
    
    cur_sac_target = rcur_sac_start_inds + info_slags(ii);
    uu = find(cur_sac_target > 1 & cur_sac_target < NT);
    cur_sac_target = cur_sac_target(uu);
    cur_sac_target(all_trialvec(used_inds(cur_sac_target)) ~= rsaccade_stop_trial_inds(uu)) = [];
    Xsac_rstart(cur_sac_target,ii) = 1;
end
Xsac_start = Xsac_start(cc_uinds,:);
Xsac_rstart = Xsac_rstart(cc_uinds,:);

cur_sac_durs = cur_sac_stop_inds - cur_sac_start_inds;


% tsince_prev_sacstart = -Inf*ones(NT,1);
% tsince_prev_fixstart = nan(NT,1);
% for ii = 1:length(cur_sac_start_inds)-1
%     cur_inds = cur_sac_start_inds(ii):cur_sac_start_inds(ii+1);
%     next_tstop = trial_end_inds(find(trial_end_inds >= cur_sac_start_inds(ii),1));
%     cur_inds(cur_inds > next_tstop) = [];
%     tsince_prev_sacstart(cur_inds) = (1:length(cur_inds))-1;
%    
%     cur_inds = cur_sac_stop_inds(ii):cur_sac_stop_inds(ii+1);
%     next_tstop = trial_end_inds(find(trial_end_inds >= cur_sac_stop_inds(ii),1));
%     cur_inds(cur_inds > next_tstop) = [];
%     tsince_prev_fixstart(cur_inds) = (1:length(cur_inds))-1;
% end
% tsince_prev_sacstart = tsince_prev_sacstart(cc_uinds);
% tsince_prev_fixstart = tsince_prev_fixstart(cc_uinds);

%initialize linear model
[~,Tinds] = meshgrid(1/spatial_usfac:1/spatial_usfac:full_nPix,1:flen);
Tinds = Tinds(use_kInds_up);

%% NOW DO INFO TIMING CALCS
%create a shuffled version of the stimulus, where different frames
%(within the original used_inds set) are randomly shuffled with
%replacement
shuf_stim = all_shift_stimmat_up;
shuf_stim(used_inds,:) = all_shift_stimmat_up(used_inds(randi(NT,NT,1)),:);
shuf_X = create_time_embedding(shuf_stim,stim_params_us);
shuf_X = shuf_X(used_inds(cc_uinds),use_kInds_up);

sacGainMod = sacStimProc(cc).gsacGainMod;
[gainLL,gain_pred_rate,G] = eval_sacgain_mod( sacGainMod, cur_Robs, all_Xmat_shift, cur_Xsac);
gain_LLseq = cur_Robs.*log2(gain_pred_rate)-gain_pred_rate;

[base_LL,~,base_predrate] = NMMmodel_eval(cur_GQM,cur_Robs,all_Xmat_shift);
base_LLseq = cur_Robs.*log2(base_predrate)-base_predrate;

temp_sparams = NMMcreate_stim_params(1);
tempmod = NMMinitialize_model(temp_sparams,1,{'lin'});

% is_before = bsxfun(@minus,tsince_prev_sacstart,Tinds') < 0;
% is_after = bsxfun(@minus,tsince_prev_fixstart,Tinds') >= 0;


% identify indices occuring within saccades
all_during = [];
for ss = 1:length(cur_sac_start_inds)
    cur_during = (cur_sac_start_inds(ss)):(cur_sac_stop_inds(ss)-1);
    all_during = cat(2,all_during,cur_during);
end
is_during = false(length(all_t_axis),1);
is_during(used_inds(all_during)) = true;

%create time-embedded Xmat. Logical indicator of isduring sac
Xis_during = repmat(is_during,[1 use_nPix_us]);
Xis_during = create_time_embedding(Xis_during,NMMcreate_stim_params([flen use_nPix_us]));
Xis_during = logical(Xis_during(used_inds(cc_uinds),:));

all_during = [];
for ss = 1:length(rcur_sac_start_inds)
    cur_during = (rcur_sac_start_inds(ss)):(rcur_sac_stop_inds(ss)-1);
    all_during = cat(2,all_during,cur_during);
end
ris_during = false(length(all_t_axis),1);
ris_during(used_inds(all_during)) = true;

rXis_during = repmat(ris_during,[1 use_nPix_us]);
rXis_during = create_time_embedding(rXis_during,NMMcreate_stim_params([flen use_nPix_us]));
rXis_during = logical(rXis_during(used_inds(cc_uinds),:));

%%
for ii = 1:length(info_slags)
    cur_set = find(Xsac_end(:,ii) == 1);
    cur_X = all_Xmat_shift(cur_set,:);
    
    [unshuf_LL,unshuf_pred_rate,G] = eval_sacgain_mod( sacGainMod, cur_Robs(cur_set), cur_X, cur_Xsac(cur_set,:));
    tempmod = NMMfit_filters(tempmod,cur_Robs(cur_set),G);
    [~,~,tempprate] = NMMmodel_eval(tempmod,cur_Robs(cur_set),G);
    unshuf_LLseq = cur_Robs(cur_set).*log2(tempprate)-tempprate;
    
    scramb_lags = find(Tinds-1 > info_slags(ii));
    cur_X(:,scramb_lags) = shuf_X(cur_set,scramb_lags);
    cur_X(Xis_during(cur_set,:)) = shuf_X(Xis_during(cur_set,:));
    
    [shuf_LL,shuf_pred_rate,Gshuff] = eval_sacgain_mod( sacGainMod, cur_Robs(cur_set), cur_X, cur_Xsac(cur_set,:));
    tempmod = NMMfit_filters(tempmod,cur_Robs(cur_set),Gshuff);
    [~,~,tempprate] = NMMmodel_eval(tempmod,cur_Robs(cur_set),Gshuff);
    shuf_LLseq = cur_Robs(cur_set).*log2(tempprate)-tempprate;
    
    %         LL_imp_before(rr,ii) = nansum((gain_LLseq(cur_set) - shuf_LLseq))/sum(cur_Robs(cur_set));
    LL_imp_before(rr,ii) = nansum(unshuf_LLseq - shuf_LLseq)/sum(cur_Robs(cur_set));
    
    cur_set = find(Xsac_rend(:,ii) == 1);
    cur_X = all_Xmat_shift(cur_set,:);
    
    [base_shufLL,~,base_shuff_prate,~,~,temp_fgint] = NMMmodel_eval(cur_GQM,cur_Robs(cur_set),cur_X);
    G = sum(bsxfun(@times,temp_fgint,[cur_GQM.mods(:).sign]),2);
    tempmod = NMMfit_filters(tempmod,cur_Robs(cur_set),G);
    [~,~,tempprate] = NMMmodel_eval(tempmod,cur_Robs(cur_set),G);
    base_unshuf_LLseq = cur_Robs(cur_set).*log2(tempprate)-tempprate;
    
    scramb_lags = find(Tinds-1 > info_slags(ii));
    cur_X(:,scramb_lags) = shuf_X(cur_set,scramb_lags);
    cur_X(rXis_during(cur_set,:)) = shuf_X(rXis_during(cur_set,:));
    
    [base_shufLL,~,base_shuff_prate,~,~,temp_fgint] = NMMmodel_eval(cur_GQM,cur_Robs(cur_set),cur_X);
    Gshuff = sum(bsxfun(@times,temp_fgint,[cur_GQM.mods(:).sign]),2);
    tempmod = NMMfit_filters(tempmod,cur_Robs(cur_set),Gshuff);
    [~,~,tempprate] = NMMmodel_eval(tempmod,cur_Robs(cur_set),Gshuff);
    base_shuf_LLseq = cur_Robs(cur_set).*log2(tempprate)-tempprate;
    
    %         base_LL_imp_before(rr,ii) = nansum((base_LLseq(cur_set) - base_shuf_LLseq))/sum(cur_Robs(cur_set));
    base_LL_imp_before(rr,ii) = nansum(base_unshuf_LLseq - base_shuf_LLseq)/sum(cur_Robs(cur_set));
end

%%
for ii = 1:length(info_slags)
    cur_set = find(Xsac_end(:,ii) == 1);
    cur_X = all_Xmat_shift(cur_set,:);
    
    [shuf_LL,shuf_pred_rate,G] = eval_sacgain_mod( sacGainMod, cur_Robs(cur_set), cur_X, cur_Xsac(cur_set,:));
    tempmod = NMMfit_filters(tempmod,cur_Robs(cur_set),G);
    [~,~,tempprate] = NMMmodel_eval(tempmod,cur_Robs(cur_set),G);
    unshuf_LLseq = cur_Robs(cur_set).*log2(tempprate)-tempprate;
    
    scramb_lags = find((Tinds-1) <= info_slags(ii));
    cur_X(:,scramb_lags) = shuf_X(cur_set,scramb_lags);
    
    [shuf_LL,shuf_pred_rate,Gshuff] = eval_sacgain_mod( sacGainMod, cur_Robs(cur_set), cur_X, cur_Xsac(cur_set,:));
    tempmod = NMMfit_filters(tempmod,cur_Robs(cur_set),Gshuff);
    [~,~,tempprate] = NMMmodel_eval(tempmod,cur_Robs(cur_set),Gshuff);
    shuf_LLseq = cur_Robs(cur_set).*log2(tempprate)-tempprate;
    
    %         LL_imp_after(rr,ii) = nansum((gain_LLseq(cur_set) - shuf_LLseq))/sum(cur_Robs(cur_set));
    LL_imp_after(rr,ii) = nansum(unshuf_LLseq - shuf_LLseq)/sum(cur_Robs(cur_set));
    
    
    cur_set = find(Xsac_rend(:,ii) == 1);
    cur_X = all_Xmat_shift(cur_set,:);
    
    [base_LL,~,base_predrate] = NMMmodel_eval(cur_GQM,cur_Robs(cur_set),all_Xmat_shift(cur_set,:));
    [base_shufLL,~,base_shuff_prate,~,~,temp_fgint] = NMMmodel_eval(cur_GQM,cur_Robs(cur_set),cur_X);
    G = sum(bsxfun(@times,temp_fgint,[cur_GQM.mods(:).sign]),2);
    tempmod = NMMfit_filters(tempmod,cur_Robs(cur_set),G);
    [~,~,tempprate] = NMMmodel_eval(tempmod,cur_Robs(cur_set),G);
    base_unshuf_LLseq = cur_Robs(cur_set).*log2(tempprate)-tempprate;
    
    scramb_lags = find((Tinds-1) <= info_slags(ii));
    cur_X(:,scramb_lags) = shuf_X(cur_set,scramb_lags);
    
    [base_LL,~,base_predrate] = NMMmodel_eval(cur_GQM,cur_Robs(cur_set),all_Xmat_shift(cur_set,:));
    [base_shufLL,~,base_shuff_prate,~,~,temp_fgint] = NMMmodel_eval(cur_GQM,cur_Robs(cur_set),cur_X);
    Gshuff = sum(bsxfun(@times,temp_fgint,[cur_GQM.mods(:).sign]),2);
    tempmod = NMMfit_filters(tempmod,cur_Robs(cur_set),Gshuff);
    [~,~,tempprate] = NMMmodel_eval(tempmod,cur_Robs(cur_set),Gshuff);
    base_shuf_LLseq = cur_Robs(cur_set).*log2(tempprate)-tempprate;
    
    %         base_LL_imp_after(rr,ii) = nansum((base_LLseq(cur_set) - base_shuf_LLseq))/sum(cur_Robs(cur_set));
    base_LL_imp_after(rr,ii) = nansum(base_unshuf_LLseq - base_shuf_LLseq)/sum(cur_Robs(cur_set));
    
end


%%
cur_X = all_Xmat_shift;
cur_X(Xis_during) = shuf_X(Xis_during);
cur_rX = all_Xmat_shift;
cur_rX(rXis_during) = shuf_X(rXis_during);
%
for ii = 1:length(info_slags)
    cur_set = find(Xsac_end(:,ii) == 1);
    
    [shuf_LL,shuf_pred_rate,G] = eval_sacgain_mod( sacGainMod, cur_Robs(cur_set), all_Xmat_shift(cur_set,:), cur_Xsac(cur_set,:));
    tempmod = NMMfit_filters(tempmod,cur_Robs(cur_set),G);
    [~,~,tempprate] = NMMmodel_eval(tempmod,cur_Robs(cur_set),G);
    unshuf_LLseq = cur_Robs(cur_set).*log2(tempprate)-tempprate;
    
    [shuf_LL,shuf_pred_rate,Gshuff] = eval_sacgain_mod( sacGainMod, cur_Robs(cur_set), cur_X(cur_set,:), cur_Xsac(cur_set,:));
    tempmod = NMMfit_filters(tempmod,cur_Robs(cur_set),Gshuff);
    [~,~,tempprate] = NMMmodel_eval(tempmod,cur_Robs(cur_set),Gshuff);
    shuf_LLseq = cur_Robs(cur_set).*log2(tempprate)-tempprate;
    
    %         LL_imp_during(rr,ii) = nansum((gain_LLseq(cur_set) - shuf_LLseq))/sum(cur_Robs(cur_set));
    LL_imp_during(rr,ii) = nansum(unshuf_LLseq - shuf_LLseq)/sum(cur_Robs(cur_set));
    
    cur_set = find(Xsac_rend(:,ii) == 1);
    
    [base_shufLL,~,base_shuff_prate,~,~,temp_fgint] = NMMmodel_eval(cur_GQM,cur_Robs(cur_set),all_Xmat_shift(cur_set,:));
    G = sum(bsxfun(@times,temp_fgint,[cur_GQM.mods(:).sign]),2);
    tempmod = NMMfit_filters(tempmod,cur_Robs(cur_set),G);
    [~,~,tempprate] = NMMmodel_eval(tempmod,cur_Robs(cur_set),G);
    base_unshuf_LLseq = cur_Robs(cur_set).*log2(tempprate)-tempprate;
    
    [base_shufLL,~,base_shuff_prate,~,~,temp_fgint] = NMMmodel_eval(cur_GQM,cur_Robs(cur_set),cur_rX(cur_set,:));
    Gshuff = sum(bsxfun(@times,temp_fgint,[cur_GQM.mods(:).sign]),2);
    tempmod = NMMfit_filters(tempmod,cur_Robs(cur_set),Gshuff);
    [~,~,tempprate] = NMMmodel_eval(tempmod,cur_Robs(cur_set),Gshuff);
    base_shuf_LLseq = cur_Robs(cur_set).*log2(tempprate)-tempprate;
    
    %         base_LL_imp_during(rr,ii) = nansum((base_LLseq(cur_set) - base_shuf_LLseq))/sum(cur_Robs(cur_set));
    base_LL_imp_during(rr,ii) = nansum(base_unshuf_LLseq - base_shuf_LLseq)/sum(cur_Robs(cur_set));
    
end
