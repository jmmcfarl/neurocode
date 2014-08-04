%% NOW INFER DRIFT CORRECTIONS
fin_fix_corr = nan(NT,1);
fin_fix_std = nan(NT,1);
fin_fix_corr(~isnan(fix_ids)) = it_fix_post_mean(end,fix_ids(~isnan(fix_ids)));
fin_fix_corr = interp1(find(~isnan(fix_ids)),fin_fix_corr(~isnan(fix_ids)),1:NT);
fin_fix_std(~isnan(fix_ids)) = it_fix_post_std(end,fix_ids(~isnan(fix_ids)));
fin_fix_std = interp1(find(~isnan(fix_ids)),fin_fix_std(~isnan(fix_ids)),1:NT);

fin_fix_corr = fin_fix_corr*sp_dx;
fin_fix_std = fin_fix_std*sp_dx;

usable_inds = find(~ismember(all_blockvec(used_inds),ignore_blocks));

measured_seqL = corrected_eye_vals_interp(used_inds,2);
measured_seqR = corrected_eye_vals_interp(used_inds,4);

min_fix_dur = 0.15;
measured_driftL = nan(size(fin_tot_corr));
measured_driftR = nan(size(fin_tot_corr));
measured_fix_avgL = nan(n_fixs,1);
measured_fix_avgR = nan(n_fixs,1);
for ii = 1:n_fixs
    cur_inds = fix_start_inds(ii):fix_stop_inds(ii);
    if length(cur_inds)*dt >= min_fix_dur & ismember(used_inds(cur_inds(1)),usable_inds)
        measured_fix_avgL(ii) = median(measured_seqL(cur_inds));
        measured_fix_avgR(ii) = median(measured_seqR(cur_inds));
        measured_driftL(cur_inds) = measured_seqL(cur_inds) - measured_fix_avgL(ii);
        measured_driftR(cur_inds) = measured_seqR(cur_inds) - measured_fix_avgR(ii);
    end
end

%overall prior on shifts
leps_prior = -(shifts*sp_dx).^2/(2*drift_jump_sigma^2);
leps_prior = bsxfun(@minus,leps_prior,logsumexp(leps_prior)); %normalize
lA_tflip = repmat(leps_prior,n_shifts,1);

usu_inds = su_inds([1 3 4 7 8]);

poss_drift_sigmas = [0.0075 0.01 0.0125 0.015 0.02];
for xv = 1:length(usu_inds)
    fprintf('Inferring drift corrections, XV %d of %d\n',xv,length(usu_inds));
    
    %% PREPROCESS MODEL COMPONENTS
    cur_xv_cell = tr_set(usu_inds(xv));
    cur_tr_cells = setdiff(tr_set,cur_xv_cell);
    cur_tr_cinds = find(ismember(tr_set,cur_tr_cells));
    cur_xv_cind = find(ismember(tr_set,cur_xv_cell));
    cur_n_tr_chs = length(cur_tr_cells);
    
    filt_bank = zeros(n_squared_filts+1,klen_us);
    cur_Xtargs = [dit_mods{nn}(cur_xv_cell).mods(:).Xtarget];
    cur_k = [dit_mods{nn}(cur_xv_cell).mods(cur_Xtargs == 1).filtK];
    n_used_filts = size(cur_k,2);
    filt_bank(1:n_used_filts,:) = cur_k';
    mod_spkNL_params = dit_mods_spkNL{nn}(cur_xv_cell).spk_NL_params;
    lin_kerns = dit_mods{nn}(cur_xv_cell).mods(cur_Xtargs == 2).filtK;
    if use_sac_kerns
        sac_kerns = dit_mods{nn}(cur_xv_cell).mods(cur_Xtargs == 3).filtK;
        msac_kerns = dit_mods{nn}(cur_xv_cell).mods(cur_Xtargs == 4).filtK;
    end
    
    %indicator predictions
    block_out = Xblock(used_inds,:)*lin_kerns;
    if use_sac_kerns
        sac_out = Xsac*sac_kerns;
        msac_out = Xmsac*msac_kerns;
    end
    %% ESTIMATE LL for each shift in each stimulus frame
    
    %precompute LL at all shifts for all units
    xv_frame_LLs = nan(NT,n_shifts);
    for xx = 1:length(shifts)
        fprintf('Shift %d of %d\n',xx,n_shifts);
        cur_stim_shift = all_Xmat_up_fixcor*shift_mat{xx};
        
        %outputs of stimulus models at current X-matrix shift
        gfuns = ones(NT,1)*mod_spkNL_params(1);
        gfuns = gfuns + cur_stim_shift*filt_bank(1,:)';
        for ff = 2:(n_squared_filts+1)
            gfuns = gfuns + mod_signs(ff)*(cur_stim_shift*filt_bank(ff,:)').^2;
        end
        
        %add contributions from extra lin kernels
        gfuns = gfuns + block_out;
        if use_sac_kerns
            gfuns = gfuns + sac_out + msac_out;
        end
        
        %incorporate beta
        gfuns = gfuns*mod_spkNL_params(2);
        
        %handle numerical overflow with log(1+exp)
        too_large = gfuns > 50;
        pred_rate = log(1+exp(gfuns));
        pred_rate(too_large) = gfuns(too_large);
        
        %incorporate alpha
        pred_rate = pred_rate*mod_spkNL_params(3);
        
        %enforce min predicted rate
        pred_rate(pred_rate < 1e-50) = 1e-50;
        
        xv_frame_LLs(:,xx) = squeeze(nansum(Robs_mat(:,cur_xv_cind).*log(pred_rate) - pred_rate,2));
    end
    
    cur_frame_LLs = frame_LLs - xv_frame_LLs; %subtract LL contribution of current LOOXV cell
    
    %% INFER MICRO-SAC SEQUENCE
    
    for pp = 1:length(poss_drift_sigmas)
        drift_sigma = poss_drift_sigmas(pp);
        
        cdist = squareform(pdist(shifts'*sp_dx));
        lA = -cdist.^2/(2*drift_sigma^2);
        lA = bsxfun(@minus,lA,logsumexp(lA,2)); %normalize
                
        lalpha=zeros(NT,n_shifts);
        lbeta = zeros(NT,n_shifts);
        %compute rescaled forward messages
        lalpha(1,:) = leps_prior + cur_frame_LLs(1,:);
        for t=2:NT
            if use_prior(t)
                cur_lA = lA_tflip;
            else
                cur_lA = lA;
            end
            lalpha(t,:) = logmulexp(lalpha(t-1,:),cur_lA) + cur_frame_LLs(t,:);
        end
        
        %compute rescaled backward messages
        lbeta(NT,:)=log(ones(1,n_shifts));
        for t=NT-1:-1:1
            if use_prior(t+1)
                cur_lA = lA_tflip;
            else
                cur_lA = lA;
            end
            lf1 = lbeta(t+1,:) + cur_frame_LLs(t+1,:);
            lbeta(t,:) = logmulexp(lf1,cur_lA');
        end
        lgamma= lalpha + lbeta;
        lgamma = bsxfun(@minus,lgamma,logsumexp(lgamma,2));
        
        gamma = exp(lgamma);
        drift_post_mean_test(xv,pp,nn,:) = sum(bsxfun(@times,gamma,shifts),2);
        drift_post_mean_cor = squeeze(drift_post_mean_test(xv,pp,nn,:));
        
        %% RECOMPUTE XMAT
        drift_Xmat = reshape(all_Xmat_up_fixcor,[NT flen full_nPix_us]);
        for ii = 1:NT
            drift_Xmat(ii,:,:) = shift_matrix_Nd(drift_Xmat(ii,:,:),-round(drift_post_mean_cor(ii)),3);
        end
        drift_Xmat = reshape(drift_Xmat,[NT flen*full_nPix_us]);
        
        %% REFIT XV CELLS
        cur_X{1} = drift_Xmat(:,use_kInds_up);
        
        silent = 1;
        cur_cell = cur_xv_cell;
        fprintf('Refitting model \n');
        cur_unit_ind = find(tr_set == cur_cell);
        cur_tr_uset = tr_inds(~isnan(Robs_mat(tr_inds,cur_unit_ind)));
        cur_xv_uset = xv_inds(~isnan(Robs_mat(xv_inds,cur_unit_ind)));
        
        tr_X = get_Xcell_tInds(cur_X,cur_tr_uset);
        xv_X = get_Xcell_tInds(cur_X,cur_xv_uset);
        
        dit_mods_test{nn+1}(cur_cell) = dit_mods_LOO{nn}(cur_cell);
        dit_mods_test{nn+1}(cur_cell) = NMMfit_filters(dit_mods_test{nn+1}(cur_cell),Robs_mat(cur_tr_uset,cur_unit_ind),...
            tr_X,[],[],silent); %fit stimulus filters
        
        %refit spk NL
        dit_mods_spkNL_test{nn+1}(cur_cell) = NMMfit_logexp_spkNL(dit_mods_test{nn+1}(cur_cell),Robs_mat(cur_tr_uset,cur_unit_ind),tr_X);
        
        newLL = NMMmodel_eval(dit_mods_spkNL_test{nn+1}(cur_cell),Robs_mat(cur_tr_uset,cur_unit_ind),tr_X);
        newxvLL = NMMmodel_eval(dit_mods_spkNL_test{nn+1}(cur_cell),Robs_mat(cur_xv_uset,cur_unit_ind),xv_X);
        
        dit_LLimp_test(pp,cur_cell) = (newLL - null_LL(cur_cell))/log(2);
        dit_xvLLimp_test(pp,cur_cell) = (newxvLL - null_xvLL(cur_cell))/log(2);
        
        %%
        
        fin_drift_corr = drift_post_mean(end,:)*sp_dx;
        fin_drift_std = drift_post_std(end,:)*sp_dx;
        for ii = 1:n_fixs
            cur_inds = fix_start_inds(ii):fix_stop_inds(ii);
            if length(cur_inds) > sac_shift
                fin_drift_corr(cur_inds(1:end-sac_shift+1)) = fin_drift_corr(cur_inds(sac_shift:end));
                fin_drift_std(cur_inds(1:end-sac_shift+1)) = fin_drift_std(cur_inds(sac_shift:end));
            end
        end
        
        fin_tot_corr = fin_fix_corr + fin_drift_corr;
        fin_tot_std = sqrt(fin_fix_std.^2 + fin_drift_std.^2);
        
        inferred_drift = nan(size(fin_tot_corr));
        inferred_fix_avg = nan(n_fixs,1);
        for ii = 1:n_fixs
            cur_inds = fix_start_inds(ii):fix_stop_inds(ii);
            if length(cur_inds)*dt >= min_fix_dur & ismember(used_inds(cur_inds(1)),usable_inds)
                cur_inf = fin_tot_corr(cur_inds);
                inferred_fix_avg(ii) = median(fin_tot_corr(cur_inds));
                inferred_drift(cur_inds) = cur_inf - inferred_fix_avg(ii);
            end
        end
        
        u = find(~isnan(measured_driftL) & ~isnan(inferred_drift));
        [drift_corrs,drif_pvals] = corr([measured_driftL(u)' measured_driftR(u)'],inferred_drift(u)','type','spearman');
        [drift_corrs_pear,drif_pvals] = corr([measured_driftL(u)' measured_driftR(u)'],inferred_drift(u)');
        [tot_corrs,tot_pvals] = corr([measured_seqL(usable_inds) measured_seqR(usable_inds)],fin_tot_corr(usable_inds)','type','spearman');
        [tot_corrs_pear,tot_pvals_pear] = corr([measured_seqL(usable_inds) measured_seqR(usable_inds)],fin_tot_corr(usable_inds)','type','pearson');
        
        drift_corr_test(xv,pp) = drift_corrs(1);
        drift_corr_pear_test(xv,pp) = drift_corrs_pear(1);
        tot_corrs_test(xv,pp) = tot_corrs(1);
        tot_corrs_pear_test(xv,pp) = tot_corrs_pear(1);
        
    end
    
end
