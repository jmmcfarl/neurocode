
info_backlag = round(0.3/dt);
info_forlag = round(0.3/dt);
info_slags = -info_backlag:info_forlag;

%lag index values in Xmat
[~,Tinds] = meshgrid(1/spatial_usfac:1/spatial_usfac:full_nPix,1:flen);
Tinds = Tinds(use_kInds_up);


cur_sac_start_inds = saccade_start_inds(big_sacs);
cur_sac_stop_inds = saccade_stop_inds(big_sacs);

% badsacs = find(cur_sac_start_inds(2:end) - cur_sac_stop_inds(1:end-1) < 50);
% cur_sac_start_inds(badsacs+1) = [];
% cur_sac_stop_inds(badsacs+1) = [];

saccade_stop_trial_inds = all_trialvec(used_inds(cur_sac_stop_inds));
cur_sac_durs = cur_sac_stop_inds - cur_sac_start_inds;

sac_jit_amp = 100;
rand_jit = randi(sac_jit_amp,length(cur_sac_start_inds),1);
rcur_sac_start_inds = cur_sac_start_inds + rand_jit;
rcur_sac_stop_inds = cur_sac_stop_inds + rand_jit;
bad = find(rcur_sac_stop_inds > NT);
rcur_sac_start_inds(bad) = [];
rcur_sac_stop_inds(bad) = [];
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


%store duration of previous saccade
prev_sac_dur = ones(NT,1)*sac_durs(end);
rprev_sac_dur = ones(NT,1)*sac_durs(end);
for ii = 1:length(cur_sac_start_inds)-1
    cur_inds = cur_sac_start_inds(ii):cur_sac_start_inds(ii+1);
    next_tstop = trial_end_inds(find(trial_end_inds >= cur_sac_start_inds(ii),1));
    cur_inds(cur_inds > next_tstop) = [];
    prev_sac_dur(cur_inds) = sac_durs(ii);

    cur_inds = rcur_sac_start_inds(ii):rcur_sac_start_inds(ii+1);
    next_tstop = trial_end_inds(find(trial_end_inds >= rcur_sac_start_inds(ii),1));
    cur_inds(cur_inds > next_tstop) = [];
    rprev_sac_dur(cur_inds) = sac_durs(ii);
end
prev_sac_dur = round(prev_sac_dur(cc_uinds)/dt);
rprev_sac_dur = round(rprev_sac_dur(cc_uinds)/dt);


%% NOW DO INFO TIMING CALCS
%create a shuffled version of the stimulus, where different frames
%(within the original used_inds set) are randomly shuffled with
%replacement
shuf_stim = all_shift_stimmat_up;
shuf_stim(used_inds,:) = all_shift_stimmat_up(used_inds(randi(NT,NT,1)),:);
shuf_X = create_time_embedding(shuf_stim,stim_params_us);
shuf_X = shuf_X(used_inds(cc_uinds),use_kInds_up);

sacGainMod = sacStimProc(cc).gsacGainMod;
% [gainLL,gain_pred_rate,G] = eval_sacgain_mod( sacGainMod, cur_Robs, all_Xmat_shift, cur_Xsac);
% gain_LLseq = cur_Robs.*log2(gain_pred_rate)-gain_pred_rate;
% 
% [base_LL,~,base_predrate] = NMMmodel_eval(cur_GQM,cur_Robs,all_Xmat_shift);
% base_LLseq = cur_Robs.*log2(base_predrate)-base_predrate;

%initialize simple 1-parameter GLM for gain/offset fitting
temp_sparams = NMMcreate_stim_params(1);
tempmod = NMMinitialize_model(temp_sparams,1,{'lin'});

%%
for ii = 1:length(info_slags)
    cur_set = find(Xsac_end(:,ii) == 1);
    cur_X = all_Xmat_shift(cur_set,:);
    
    [unshuf_LL,unshuf_pred_rate,G] = eval_sacgain_mod( sacGainMod, cur_Robs(cur_set), cur_X, cur_Xsac(cur_set,:));
    tempmod = NMMfit_filters(tempmod,cur_Robs(cur_set),G);
    [~,~,tempprate] = NMMmodel_eval(tempmod,cur_Robs(cur_set),G);
    unshuf_LLseq = cur_Robs(cur_set).*log2(tempprate)-tempprate;
    
    
    
    cur_Tinds = bsxfun(@minus,prev_sac_dur(cur_set),Tinds');
    scramb_lags = find(cur_Tinds-1 >= info_slags(ii));
    cur_shufX = shuf_X(cur_set,:);
    scr_X = cur_X;
    scr_X(scramb_lags) = cur_shufX(scramb_lags);
    [shuf_LL,shuf_pred_rate,Gshuff] = eval_sacgain_mod( sacGainMod, cur_Robs(cur_set), scr_X, cur_Xsac(cur_set,:));
    tempmod = NMMfit_filters(tempmod,cur_Robs(cur_set),Gshuff);
    [~,~,tempprate] = NMMmodel_eval(tempmod,cur_Robs(cur_set),Gshuff);
    shuf_LLseq = cur_Robs(cur_set).*log2(tempprate)-tempprate;
    LL_imp_before(ii) = nansum(unshuf_LLseq - shuf_LLseq)/sum(cur_Robs(cur_set));
    
    
    scramb_lags = find((Tinds-1) <= info_slags(ii));
    scr_X = cur_X;
    scr_X(:,scramb_lags) = shuf_X(cur_set,scramb_lags);
    [shuf_LL,shuf_pred_rate,Gshuff] = eval_sacgain_mod( sacGainMod, cur_Robs(cur_set), scr_X, cur_Xsac(cur_set,:));
    tempmod = NMMfit_filters(tempmod,cur_Robs(cur_set),Gshuff);
    [~,~,tempprate] = NMMmodel_eval(tempmod,cur_Robs(cur_set),Gshuff);
    shuf_LLseq = cur_Robs(cur_set).*log2(tempprate)-tempprate;
    LL_imp_after(ii) = nansum(unshuf_LLseq - shuf_LLseq)/sum(cur_Robs(cur_set));
    
    
    poss_scramb_lags = find((Tinds-1) > info_slags(ii)); %set of lags before fixation onset
    
    %set of lags after saccade onset
    cur_Tinds = bsxfun(@minus,prev_sac_dur(cur_set),Tinds');
    scramb_lags = find(cur_Tinds-1 < info_slags(ii));
    
    [scramb_i,scramb_j] = ind2sub([length(cur_set) length(Tinds)],scramb_lags);
    scramb_lags = scramb_lags(ismember(scramb_j,poss_scramb_lags));
    
     scr_X = cur_X;
   cur_shufX = shuf_X(cur_set,:);
    scr_X(scramb_lags) = cur_shufX(scramb_lags);
    [shuf_LL,shuf_pred_rate,Gshuff] = eval_sacgain_mod( sacGainMod, cur_Robs(cur_set), scr_X, cur_Xsac(cur_set,:));
    tempmod = NMMfit_filters(tempmod,cur_Robs(cur_set),Gshuff);
    [~,~,tempprate] = NMMmodel_eval(tempmod,cur_Robs(cur_set),Gshuff);
    shuf_LLseq = cur_Robs(cur_set).*log2(tempprate)-tempprate;
    LL_imp_during(ii) = nansum(unshuf_LLseq - shuf_LLseq)/sum(cur_Robs(cur_set));

end


%%
Nrpts = 5;
sac_jit_amp = 100;
for rr = 1:Nrpts
    fprintf('Scrambling with rand saclocs, iter %d of %d\n',rr,Nrpts);
    
    rand_jit = randi(sac_jit_amp,length(cur_sac_start_inds),1);
    rcur_sac_start_inds = cur_sac_start_inds + rand_jit;
    rcur_sac_stop_inds = cur_sac_stop_inds + rand_jit;
    bad = find(rcur_sac_stop_inds > NT);
    rcur_sac_start_inds(bad) = [];
    rcur_sac_stop_inds(bad) = [];
    rsaccade_stop_trial_inds = all_trialvec(used_inds(rcur_sac_stop_inds));
    
    Xsac_rend = zeros(NT,length(info_slags));
    for ii = 1:length(info_slags)
        cur_sac_target = rcur_sac_stop_inds + info_slags(ii);
        uu = find(cur_sac_target > 1 & cur_sac_target < NT);
        cur_sac_target = cur_sac_target(uu);
        cur_sac_target(all_trialvec(used_inds(cur_sac_target)) ~= rsaccade_stop_trial_inds(uu)) = [];
        Xsac_rend(cur_sac_target,ii) = 1;
    end
    Xsac_rend = Xsac_rend(cc_uinds,:);
    
    %store duration of previous saccade
    rprev_sac_dur = ones(NT,1)*sac_durs(end);
    for ii = 1:length(cur_sac_start_inds)-1
        cur_inds = rcur_sac_start_inds(ii):rcur_sac_start_inds(ii+1);
        next_tstop = trial_end_inds(find(trial_end_inds >= rcur_sac_start_inds(ii),1));
        cur_inds(cur_inds > next_tstop) = [];
        rprev_sac_dur(cur_inds) = sac_durs(ii);
    end
    rprev_sac_dur = round(rprev_sac_dur(cc_uinds)/dt);
    
    
    
    for ii = 1:length(info_slags)
        
        cur_set = find(Xsac_rend(:,ii) == 1);
        cur_X = all_Xmat_shift(cur_set,:);
        
        [base_shufLL,~,base_shuff_prate,~,~,temp_fgint] = NMMmodel_eval(cur_GQM,cur_Robs(cur_set),cur_X);
        G = sum(bsxfun(@times,temp_fgint,[cur_GQM.mods(:).sign]),2);
        tempmod = NMMfit_filters(tempmod,cur_Robs(cur_set),G);
        [~,~,tempprate] = NMMmodel_eval(tempmod,cur_Robs(cur_set),G);
        base_unshuf_LLseq = cur_Robs(cur_set).*log2(tempprate)-tempprate;
        
        
        cur_Tinds = bsxfun(@minus,rprev_sac_dur(cur_set),Tinds');
        scramb_lags = find(cur_Tinds-1 >= info_slags(ii));
        scr_X = cur_X;
        cur_shufX = shuf_X(cur_set,:);
        scr_X(scramb_lags) = cur_shufX(scramb_lags);
        
        [base_shufLL,~,base_shuff_prate,~,~,temp_fgint] = NMMmodel_eval(cur_GQM,cur_Robs(cur_set),scr_X);
        Gshuff = sum(bsxfun(@times,temp_fgint,[cur_GQM.mods(:).sign]),2);
        tempmod = NMMfit_filters(tempmod,cur_Robs(cur_set),Gshuff);
        [~,~,tempprate] = NMMmodel_eval(tempmod,cur_Robs(cur_set),Gshuff);
        base_shuf_LLseq = cur_Robs(cur_set).*log2(tempprate)-tempprate;
        base_LL_imp_before(rr,ii) = nansum(base_unshuf_LLseq - base_shuf_LLseq)/sum(cur_Robs(cur_set));
        
        
        
        scramb_lags = find((Tinds-1) <= info_slags(ii));
        scr_X = cur_X;
        scr_X(:,scramb_lags) = shuf_X(cur_set,scramb_lags);
        
        [base_LL,~,base_predrate] = NMMmodel_eval(cur_GQM,cur_Robs(cur_set),all_Xmat_shift(cur_set,:));
        [base_shufLL,~,base_shuff_prate,~,~,temp_fgint] = NMMmodel_eval(cur_GQM,cur_Robs(cur_set),scr_X);
        Gshuff = sum(bsxfun(@times,temp_fgint,[cur_GQM.mods(:).sign]),2);
        tempmod = NMMfit_filters(tempmod,cur_Robs(cur_set),Gshuff);
        [~,~,tempprate] = NMMmodel_eval(tempmod,cur_Robs(cur_set),Gshuff);
        base_shuf_LLseq = cur_Robs(cur_set).*log2(tempprate)-tempprate;
        base_LL_imp_after(rr,ii) = nansum(base_unshuf_LLseq - base_shuf_LLseq)/sum(cur_Robs(cur_set));
        
        
        
        poss_scramb_lags = find((Tinds-1) > info_slags(ii));
        cur_Tinds = bsxfun(@minus,rprev_sac_dur(cur_set),Tinds');
        scramb_lags = find(cur_Tinds-1 < info_slags(ii));
        [scramb_i,scramb_j] = ind2sub([length(cur_set) length(Tinds)],scramb_lags);
        scramb_lags = scramb_lags(ismember(scramb_j,poss_scramb_lags));
        scr_X = cur_X;
        cur_shufX = shuf_X(cur_set,:);
        scr_X(scramb_lags) = cur_shufX(scramb_lags);
        
        [base_shufLL,~,base_shuff_prate,~,~,temp_fgint] = NMMmodel_eval(cur_GQM,cur_Robs(cur_set),scr_X);
        Gshuff = sum(bsxfun(@times,temp_fgint,[cur_GQM.mods(:).sign]),2);
        tempmod = NMMfit_filters(tempmod,cur_Robs(cur_set),Gshuff);
        [~,~,tempprate] = NMMmodel_eval(tempmod,cur_Robs(cur_set),Gshuff);
        base_shuf_LLseq = cur_Robs(cur_set).*log2(tempprate)-tempprate;
        base_LL_imp_during(rr,ii) = nansum(base_unshuf_LLseq - base_shuf_LLseq)/sum(cur_Robs(cur_set));
        
    end
end
