
%%
close all
%log prior dist on drift-corrections
drift_jump_Lprior = -(drift_shifts*drift_stim_dx).^2./(2*HMM_params.drift_jump_sigma^2); %log-prior on initial position (after a jump) during drift-inference
drift_jump_Lprior = drift_jump_Lprior - logsumexp(drift_jump_Lprior);

drift_prior_sigma = 0.0075;
%log prior matrix for drift (not including coil-signals)
cdist = squareform(pdist(drift_shifts'*drift_stim_dx));
base_LA = -cdist.^2/(2*drift_prior_sigma^2); %baseline log-prior matrix on drift position changes
base_LA = bsxfun(@minus,base_LA,logsumexp(base_LA,2)); %normalize

acc_prior_sigma = 0.01;
[XX,YY,ZZ] = ndgrid(drift_shifts*drift_stim_dx);
sec_deriv = XX - 2*YY + ZZ;
LA = -sec_deriv.^2/(2*acc_prior_sigma^2);
LA = bsxfun(@minus,LA,logsumexp(LA,3)); %normalize

base_prior = bsxfun(@plus,drift_jump_Lprior',drift_jump_Lprior);


vdrift_mean = nan(NT,1);
vdrift_std = nan(NT,1);
adrift_mean = nan(NT,1);
adrift_std = nan(NT,1);
for ff = 1:1000;
    ff
    % ff = 30
    tset = find(pfix_ids==ff)'; %indices in this projected fixation
    ntset = length(tset);
    if ntset > 1
        [talpha,tbeta] = deal(zeros(ntset,n_drift_shifts)); %initialize forward and backward messages
        cur_LL_set = frame_LLs(tset,:);
        
        %forward messages
        talpha(1,:) = drift_jump_Lprior + cur_LL_set(1,:); %initialize forward messages within this fixation
        for t = 2:ntset
            talpha(t,:) = logmulexp(talpha(t-1,:),base_LA) + cur_LL_set(t,:);
        end
        %     for_probs = bsxfun(@minus,talpha,logsumexp(talpha,2));
        %     for_probs = exp(for_probs);
        
        %backward messages
        tbeta(end,:) = zeros(1,n_drift_shifts);
        for t = (ntset-1):-1:1
            lf1 = tbeta(t+1,:) + cur_LL_set(t+1,:);
            tbeta(t,:) = logmulexp(lf1,base_LA');
        end
        temp_gamma = talpha + tbeta;
        %     back_probs = bsxfun(@minus,tbeta,logsumexp(tbeta,2));
        %     back_probs = exp(back_probs);
        
        temp_gamma = bsxfun(@minus,temp_gamma,logsumexp(temp_gamma,2));
        gamma = exp(temp_gamma);
        vdrift_mean(tset) = sum(bsxfun(@times,gamma,drift_shifts),2);
        vdrift_std(tset) = sqrt(sum(bsxfun(@times,gamma,drift_shifts.^2),2) - vdrift_mean(tset).^2);
        
        % % close all
        % f1 = figure();
        % subplot(3,1,1)
        % imagesc(for_probs');
        % subplot(3,1,2)
        % imagesc(back_probs');
        % subplot(3,1,3)
        % imagesc(gamma');
        
        %     f2 = figure; hold on
        %     shadedErrorBar(1:ntset,post_mean*drift_stim_dx,post_std*drift_stim_dx,{'color','r'});
        
        % % pause
        % % close all
        % % end
        % %
        % %initialize variables
        % psi=zeros(ntset,n_drift_shifts); %Psi stores the most likely preceeding state
        % delta = zeros(ntset,n_drift_shifts);
        %
        % [delta(1,:),psi(1,:)] = max(talpha(1,:));
        %
        % for t = 2:ntset
        %     for d = 1:n_drift_shifts
        %         temp = delta(t-1,:) + base_LA(d,:);
        %         [delta(t,d),psi(t,d)] = max(temp);
        %     end
        %     delta(t,:) = delta(t,:) + cur_LL_set(t,:);
        % end
        %
        % % Backtracking Viterbi
        % state_seq = zeros(ntset,1);
        % [llik_best,state_seq(ntset)] = max(delta(ntset,:));
        % for t=ntset-1:-1:1
        %     state_seq(t) = psi(t+1,state_seq(t+1));
        % end
        %
        % hold on
        % plot(drift_shifts(state_seq)*drift_stim_dx,'b','linewidth',2)
        %
        %
        %     drift_init_sigma = HMM_params.drift_jump_sigma;
        %     %log prior dist on drift-corrections
        %     drift_jump_Lprior = -(drift_shifts*drift_stim_dx).^2./(2*drift_init_sigma^2); %log-prior on initial position (after a jump) during drift-inference
        %     drift_jump_Lprior = drift_jump_Lprior - logsumexp(drift_jump_Lprior);
        
        
%         % %log prior matrix for drift (not including coil-signals)
%         % cdist = squareform(pdist(drift_shifts'*drift_stim_dx));
%         % prior2 = -cdist.^2/(2*(HMM_params.drift_prior_sigma*drift_dsf)^2); %baseline log-prior matrix on drift position changes
%         % prior2 = bsxfun(@minus,base_LA,logsumexp(base_LA,2)); %normalize
%         
%         % INFER DRIFT CORRECTIONS
%         % close all
%         Lgamma = nan(NT,n_drift_shifts);
%         reverseStr = '';
%         % for ff = 1:n_fixs
%         %     if mod(ff,100)==0
%         %         msg = sprintf('Inferring drift in fixation %d of %d\n',ff,n_fixs);
%         %         fprintf([reverseStr msg]);
%         %         reverseStr = repmat(sprintf('\b'),1,length(msg));
%         %     end
%         %     ff = 30;
%         tset = find(pfix_ids==ff)'; %indices in this projected fixation
%         ntset = length(tset);
%         cur_LL_set = frame_LLs(tset,:);
%         
%         [talpha,tbeta] = deal(zeros(ntset,n_drift_shifts,n_drift_shifts)); %initialize forward and backward messages
%         
%         %forward messages
%         talpha(1,:,:) = bsxfun(@plus,base_prior,cur_LL_set(1,:)); %initialize forward messages within this fixation
%         talpha(2,:,:) = bsxfun(@plus,base_prior,cur_LL_set(2,:)); %initialize forward messages within this fixation
%         for t = 3:ntset
%             temp = bsxfun(@plus,squeeze(talpha(t-1,:,:)),LA);
%             temp = bsxfun(@plus,temp,reshape(cur_LL_set(t,:),1,1,[]));
%             talpha(t,:,:) = logsumexp(temp,1);
%             
%             %         cur_mat = talpha(t,:,:);
%             %         cur_mat = cur_mat - max(cur_mat(:));
%             %         imagesc(squeeze(exp(cur_mat)))
%             %         t
%             %         pause
%         end
%         
%         %     for_probs = squeeze(logsumexp(talpha,2));
%         %     for_probs = bsxfun(@minus,for_probs,logsumexp(for_probs,2));
%         %     for_probs = exp(for_probs);
%         
%         %backward messages
%         tbeta(end,:,:) = zeros(1,n_drift_shifts,n_drift_shifts);
%         %     tbeta(end-1,:,:) = zeros(1,n_drift_shifts,n_drift_shifts);
%         %     tbeta(end,:,:) = base_prior; %initialize forward messages within this fixation
%         tbeta(end-1,:,:) = bsxfun(@plus,base_prior',cur_LL_set(end,:)); %initialize forward messages within this fixation
%         for t = (ntset-2):-1:1
%             temp = bsxfun(@plus,squeeze(tbeta(t+1,:,:)),reshape(cur_LL_set(t+1,:),1,[]));
%             %         temp = bsxfun(@plus,temp,permute(LA,[3 1 2]));
%             temp = bsxfun(@plus,reshape(temp,1,n_drift_shifts,n_drift_shifts),LA);
%             tbeta(t,:,:) = logsumexp(temp,3);
%         end
%         
%         %     back_probs = squeeze(logsumexp(tbeta,2));
%         %     back_probs = bsxfun(@minus,back_probs,logsumexp(back_probs,2));
%         %     back_probs = exp(back_probs);
%         
%         
%         temp_gamma = talpha + tbeta;
%         temp_gamma = squeeze(logsumexp(temp_gamma,2));
%         temp_gamma = bsxfun(@minus,temp_gamma,logsumexp(temp_gamma,2));
%         gamma = exp(temp_gamma);
%         adrift_mean(tset) = sum(bsxfun(@times,gamma,drift_shifts),2);
%         adrift_std(tset) = sqrt(sum(bsxfun(@times,gamma,drift_shifts.^2),2) - adrift_mean(tset).^2);
%         %     Lgamma(tset,:) = repmat(temp_gamma,ntset,1);
%         % end
%         %
%         % % close all
%         % f4 = figure();
%         % subplot(3,1,1)
%         % imagesc(for_probs');
%         % subplot(3,1,2)
%         % imagesc(back_probs');
%         % subplot(3,1,3)
%         % imagesc(gamma');
%         
%         % % f2 = figure;
%         % figure(f2); hold on
%         % shadedErrorBar(1:ntset,post_mean*drift_stim_dx,post_std*drift_stim_dx);
%         
%         
%         
%         % %initialize variables
%         % psi=zeros(ntset,n_drift_shifts,n_drift_shifts); %Psi stores the most likely preceeding state
%         % delta = zeros(ntset,n_drift_shifts,n_drift_shifts);
%         % for d1 = 1:n_drift_shifts
%         % [delta(1,d1,:),psi(1,d1,:)] = max(talpha(1,d1,:),[],3);
%         % end
%         %
%         % for t = 2:ntset
%         %     for d1 = 1:n_drift_shifts
%         %         temp = bsxfun(@plus,squeeze(delta(t-1,:,d1)),squeeze(LA(:,d1,:)));
%         %         [delta(t,d1,:),psi(t,d1,:)] = max(temp,[],2);
%         %     end
%         %     delta(t,:,:) = bsxfun(@plus,delta(t,:,:),cur_LL_set(t,:));
%         % end
%         %
%         % % Backtracking Viterbi
%         % state_seq = zeros(ntset,1);
%         % [llik_best,best_loc] = max(reshape(delta(ntset,:,:),[],1));
%         % [best_i,best_j] = ind2sub([n_drift_shifts n_drift_shifts],best_loc);
%         % state_seq(end) = best_i;
%         % [~,state_seq(end-1)] = max(delta(ntset-1,:,state_seq(end)));
%         % for t=(ntset-2):-1:1
%         %     state_seq(t) = psi(t+1,state_seq(t+1),state_seq(t+2));
%         % end
%         %
%         % hold on
%         % plot(drift_shifts(state_seq)*drift_stim_dx,'m','linewidth',2)
%         %
%         
%         % cur_inds = used_inds(tset) - neural_delay;
%         % uu = find(cur_inds >= 0);
%         % cur_EP = corrected_eye_vals_interp(cur_inds(uu),[2 4]);
%         % cur_EP = bsxfun(@minus,cur_EP,median(cur_EP));
%         % plot(uu,cur_EP,'b')
%         
%         
%         % pause
%         % close all
    end
end


%%
drift_up_ratio = drift_usfac/fix_usfac;

[vfin_tot_corr,vfin_tot_std] = construct_eye_position(drift_up_ratio*fix_it_post_mean(end,:),drift_up_ratio*fix_it_post_std(end,:),...
    vdrift_mean,vdrift_std,fix_ids,trial_start_inds,trial_end_inds,round(HMM_params.neural_delay/gen_data.dt));
[afin_tot_corr,afin_tot_std] = construct_eye_position(drift_up_ratio*fix_it_post_mean(end,:),drift_up_ratio*fix_it_post_std(end,:),...
    adrift_mean,adrift_std,fix_ids,trial_start_inds,trial_end_inds,round(HMM_params.neural_delay/gen_data.dt));

%%
close all

uNT = 2e4;
figure;hold on
% shadedErrorBar(1:uNT,afin_tot_corr(1:uNT)*drift_stim_dx,afin_tot_std(1:uNT)*drift_stim_dx,{'color','k'});
shadedErrorBar(1:uNT,vfin_tot_corr(1:uNT)*drift_stim_dx,vfin_tot_std(1:uNT)*drift_stim_dx,{'color','r'});
plot(1:uNT,corrected_eye_vals_interp(used_inds(1:uNT),[2]),'g','linewidth',2)
plot(1:uNT,corrected_eye_vals_interp(used_inds(1:uNT),[4]),'b','linewidth',2)

%%
buff = 2;
eye_drifts = nan(NT,2);
inf_drifts = nan(NT,1);
for ii = 1:ff
    cur_range = find(fix_ids == ii);
    if length(cur_range) > 20
        cur_range([1:buff]) = [];
        cur_range((end-buff+1):end) = [];
        cur_eye_samp = corrected_eye_vals_interp(used_inds(cur_range),[2 4]);
        cur_eye_samp(:,1) = jmm_smooth_1d_cor(cur_eye_samp(:,1),2);
        cur_eye_samp(:,2) = jmm_smooth_1d_cor(cur_eye_samp(:,2),2);
%         eye_drifts(cur_range(2:end),:) = diff(cur_eye_samp);
        eye_drifts(cur_range,:) = bsxfun(@minus,cur_eye_samp, nanmedian(cur_eye_samp));
        
%         inf_drifts(cur_range(2:end)) = drift_stim_dx*diff(afin_tot_corr(cur_range));
%          inf_drifts(cur_range) = drift_stim_dx*(afin_tot_corr(cur_range) - nanmedian(afin_tot_corr(cur_range)));
         inf_drifts(cur_range) = drift_stim_dx*(vfin_tot_corr(cur_range) - nanmedian(afin_tot_corr(cur_range)));
   end
end