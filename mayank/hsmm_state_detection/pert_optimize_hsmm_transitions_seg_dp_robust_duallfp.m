function [new_hsmm_state_seq] = pert_optimize_hsmm_transitions_seg_dp_robust_duallfp(new_obs,hmm,orig_state_seq,gamma,Fs_orig,Fs_new,pert_range)

%% currently works only for univariate nonstationary observations or
%% stationary multivariate observations
% min_mllik = -1;
min_mllik = hmm.min_mllik;

%these are the indices for the segments of data with UDS
UDS_samples = [];
for i = 1:hmm.Nsegs
    UDS_seg_inds{i} = hmm.UDS_segs(i,1):hmm.UDS_segs(i,2);
    curT(i) = length(UDS_seg_inds{i});
    UDS_samples = [UDS_samples hmm.UDS_segs(i,1):hmm.UDS_segs(i,2)];
    UDS_segind_new_start = round((hmm.UDS_segs(i,1)-1)*Fs_new/Fs_orig)+1;
    enddiff = length(new_obs)-UDS_segind_new_start+1;
    curT_new(i) = min(round(Fs_new/Fs_orig*(curT(i)+1)),enddiff);
    UDS_seg_inds_new{i} = UDS_segind_new_start:...
        UDS_segind_new_start+curT_new(i)-1;
    UDS_seg_inds_new{i}(UDS_seg_inds_new{i} > length(new_obs)) = [];
end

meantype = hmm.meantype;
for i = 1:hmm.Nsegs
    if size(gamma{i},1) < size(gamma{i},2)
        gamma{i} = gamma{i}';
    end
end

dsf = round(2016/Fs_orig);
dsf_new = round(2016/Fs_new);

[T_new,p] = size(new_obs);

A = hmm.A;

%compute state-dependent broadband amplitude
new_obs_d = downsample(new_obs,dsf/dsf_new);

if strcmp(meantype,'variable')
    meandiff = [];
    meandiff2 = [];
    for i = 1:hmm.Nsegs
        t_axis{i} = (1:curT(i))/Fs_orig;
        t_axis_new{i} = (1:curT_new(i))/Fs_new;
        total_dur = curT(i)/hmm.Fs;
        numWins = floor((total_dur-hmm.windowSize)/hmm.windowSlide);
        win_t{i} = (0:numWins-1)*hmm.windowSlide+hmm.windowSize/2;
        win_t{i} = [0 win_t{i} max(t_axis_new{i})];
        temp_meanfun{i} = zeros(2,numWins);
        temp_meanfun2{i} = zeros(2,numWins);
        tot_prob = zeros(2,numWins);
        for l = 1:2
            for w = 1:numWins
                begT = round((w-1)*hmm.windowSlide*hmm.Fs)+1;
                endT = begT + round(hmm.windowSize*hmm.Fs);
                temp_meanfun{i}(l,w) = sum(gamma{i}(begT:endT,l).*new_obs_d(UDS_seg_inds{i}(begT:endT),1))/sum(gamma{i}(begT:endT,l));
                temp_meanfun2{i}(l,w) = sum(gamma{i}(begT:endT,l).*new_obs_d(UDS_seg_inds{i}(begT:endT),2))/sum(gamma{i}(begT:endT,l));
                tot_prob(l,w) = sum(gamma{i}(begT:endT,l))/(hmm.windowSize*hmm.Fs);
            end
            temp_meanfun{i}(l,tot_prob(l,:) < 0.05) = nan;
             temp_meanfun2{i}(l,tot_prob(l,:) < 0.05) = nan;
       end
            meandiff = [meandiff temp_meanfun{i}(2,tot_prob(l,:) >= 0.05)-...
                temp_meanfun{i}(1,tot_prob(l,:) >= 0.05)];
            meandiff2 = [meandiff2 temp_meanfun2{i}(2,tot_prob(l,:) >= 0.05)-...
                temp_meanfun2{i}(1,tot_prob(l,:) >= 0.05)];
    end
    avgmeandiff = nanmean(meandiff);
    avgmeandiff2 = nanmean(meandiff2);
    for i = 1:hmm.Nsegs
        temp_meanfun{i}(1,isnan(temp_meanfun{i}(1,:))) = temp_meanfun{i}(2,isnan(temp_meanfun{i}(1,:))) - avgmeandiff;
        temp_meanfun{i}(2,isnan(temp_meanfun{i}(2,:))) = temp_meanfun{i}(1,isnan(temp_meanfun{i}(2,:))) + avgmeandiff;
        temp_meanfun{i} = [temp_meanfun{i}(:,1) temp_meanfun{i} temp_meanfun{i}(:,end)]';
        meanfun_est{i} = interp1(win_t{i},temp_meanfun{i},t_axis_new{i});
        temp_meanfun2{i}(1,isnan(temp_meanfun2{i}(1,:))) = temp_meanfun2{i}(2,isnan(temp_meanfun2{i}(1,:))) - avgmeandiff2;
        temp_meanfun2{i}(2,isnan(temp_meanfun2{i}(2,:))) = temp_meanfun2{i}(1,isnan(temp_meanfun2{i}(2,:))) + avgmeandiff2;
        temp_meanfun2{i} = [temp_meanfun2{i}(:,1) temp_meanfun2{i} temp_meanfun2{i}(:,end)]';
        meanfun_est2{i} = interp1(win_t{i},temp_meanfun2{i},t_axis_new{i});
        if t_axis{i}(2) ~= 0
            t_axis{i}(1) = 0;
        end
        if t_axis{i}(end-1) ~= t_axis_new{i}(end)
            t_axis{i}(end) = t_axis_new{i}(end);
        end
        gamma_inter{i} = interp1(t_axis{i},gamma{i},t_axis_new{i});
    end
    for l = 1:hmm.K
        o_sum = 0;
        g_sum = 0;
        for i = 1:hmm.Nsegs
            state_mean_est = [meanfun_est{i}(:,l) meanfun_est2{i}(:,l)];
            mdiff{i} = new_obs(UDS_seg_inds_new{i},:) - state_mean_est;
            o_sum = o_sum + (repmat(gamma_inter{i}(:,l),1,p).*mdiff{i})'*mdiff{i};
            g_sum = g_sum + sum(gamma_inter{i}(:,l));
        end
        var_est(l,:,:) = o_sum/g_sum;
        R_ch = chol(squeeze(var_est(l,:,:)));
        for i = 1:hmm.Nsegs
            state_mean_est = [meanfun_est{i}(:,l) meanfun_est2{i}(:,l)];
            mdiff{i} = new_obs(UDS_seg_inds_new{i},:) - state_mean_est;
            mdiff_R{i} = mrdivide(mdiff{i}, R_ch);
           logSqrtDetSigma = sum(log(diag(R_ch)));
           quadform = sum(mdiff_R{i}.^2, 2);
            state_llik{i}(l,:) = -0.5*quadform - logSqrtDetSigma - p*log(2*pi)/2;
        end
    end
    avg_var = squeeze(mean(var_est));
    R_ch = chol(avg_var);
    for i = 1:hmm.Nsegs
        max_state_llik = max(state_llik{i});
        bad_model = find(max_state_llik < min_mllik);
        if ~isempty(bad_model)
           for l = 1:hmm.K
               state_mean_est = [meanfun_est{i}(:,l) meanfun_est2{i}(:,l)];
               mdiff{i} = new_obs(UDS_seg_inds_new{i},:) - state_mean_est;
             mdiff_R{i} = mrdivide(mdiff{i}, R_ch);
           logSqrtDetSigma = sum(log(diag(R_ch)));
           quadform = sum(mdiff_R{i}.^2, 2);
             rob_llik(l,:) = -0.5*quadform - logSqrtDetSigma - p*log(2*pi)/2;
           end
           state_llik{i}(:,bad_model) = rob_llik(:,bad_model);
        end
        clear rob_llik
    end
end

%%
% eps = 1e-200;
% min_mllik = -1;
% 
% %compute likelihood vs time for each state
% for i = 1:hmm.Nsegs
%     cur_obs = new_obs(UDS_seg_inds_new{i},:);
%     curT = size(cur_obs,1);
%     B = zeros(curT,hmm.K);
%     for l=1:hmm.K
%         state_meanfun = [meanfun_est{i}(:,l) meanfun_est2{i}(:,l)];
%         state_var = hmm.state(l).var;
%         d = cur_obs-state_meanfun;
%         B(:,l) = mvnpdf(d,zeros(curT,hmm.p),state_var);
%     end
% 
%     B(B < eps) = eps;
% 
%     max_obs_lik = max(log(B),[],2);
% 
%     bad_model = find(max_obs_lik < min_mllik);
% 
%     if ~isempty(bad_model)
%         B_r = zeros(curT,hmm.K);
%         avg_var = 0;
%         for l = 1:hmm.K
%             avg_var = avg_var + hmm.state(l).var;
%         end
%         avg_var = avg_var/hmm.K;
%         for l=1:hmm.K
%             state_meanfun = [meanfun_est{i}(:,l) meanfun_est2{i}(:,l)];
%             d = cur_obs-state_meanfun;
%             B_r(:,l) = mvnpdf(d,zeros(curT,hmm.p),avg_var);
%         end
%         B(bad_model,:) = B_r(bad_model,:);
%     end
% end

%% locate state transitions in new time
pert_range = round(pert_range(1)*Fs_new):round(pert_range(2)*Fs_new);
for i = 1:hmm.Nsegs
    orig_state_seq{i} = orig_state_seq{i}(:);
    up_trans{i} = find(orig_state_seq{i}(1:end-1) == 1 & orig_state_seq{i}(2:end) == 2);
    down_trans{i} = find(orig_state_seq{i}(1:end-1) == 2 & orig_state_seq{i}(2:end) == 1);
    state_trans{i} = [up_trans{i} ones(size(up_trans{i}))*2; down_trans{i} ones(size(down_trans{i}))];
    [dummy,sort_order] = sort(state_trans{i}(:,1));
    state_trans{i} = state_trans{i}(sort_order,:);
    state_trans{i}(:,1) = round(state_trans{i}(:,1)*Fs_new/Fs_orig);
end

%%
us_fac = Fs_new/hmm.Fs;
poss_durs = linspace(hmm.min_state_dur,hmm.max_state_dur,length(hmm.dur_range)*us_fac);
for k = 1:hmm.K
    if strcmp(hmm.state(k).dur_type,'gamma')
        P(k,:) = gamma_pmf(poss_durs/Fs_new,hmm.state(k).dur_pars(1),hmm.state(k).dur_pars(2));
    elseif strcmp(hmm.state(k).dur_type,'inv_gauss')
        P(k,:) = inverse_gaussian_pmf(poss_durs/Fs_new,hmm.state(k).dur_pars(1),hmm.state(k).dur_pars(2));
    else
        P(k,:) = geometric_pmf(poss_durs,1-hmm.A(k,k));         
    end
end
%%
p_d = log(P');
num_perts = length(pert_range);
[i_mat,j_mat] = meshgrid(pert_range,pert_range);
delta_mat = i_mat - j_mat;

for i = 1:hmm.Nsegs
    psi = cumsum(state_llik{i}');
    num_trans = size(state_trans{i},1);
    length_to = -Inf*ones(num_trans,num_perts);
    coming_from = zeros(num_trans,num_perts);
    
    %for t = 1
    trans_ind = state_trans{i}(1,1);
    trans_to = state_trans{i}(1,2); 
    trans_from = setdiff([1 2],trans_to);
    poss_pert_to_inds = trans_ind + pert_range;
    poss_pert_to_range = 1:num_perts;
    poss_pert_to_range(poss_pert_to_inds <= 1 | poss_pert_to_inds > length(UDS_seg_inds_new{i})) = [];
    cur_dur = trans_ind;
    lobs_prob_pre = psi(poss_pert_to_inds(poss_pert_to_range)-1,trans_from) - psi(poss_pert_to_inds(poss_pert_to_range(1))-1,trans_from);
    lobs_prob_post = psi(poss_pert_to_inds(poss_pert_to_range(end)),trans_to) - psi(poss_pert_to_inds(poss_pert_to_range)-1,trans_to);
    lobs_prob = lobs_prob_pre + lobs_prob_post;
    prop_durs = poss_pert_to_inds(poss_pert_to_range);
    log_duration_prob = p_d(prop_durs,trans_from);
    length_to(1,poss_pert_to_range) = lobs_prob + log_duration_prob;
    
    for t = 2:num_trans     
        trans_ind = state_trans{i}(t,1);
        trans_to = state_trans{i}(t,2);
        trans_from = setdiff([1 2],trans_to);
        poss_pert_to_inds = trans_ind + pert_range;
        poss_pert_to_range = 1:num_perts;            
        poss_pert_to_range(poss_pert_to_inds <= 1 | poss_pert_to_inds > length(UDS_seg_inds_new{i})) = [];
        cur_dur = trans_ind - state_trans{i}(t-1,1);
        lobs_prob_pre = psi(poss_pert_to_inds(poss_pert_to_range)-1,trans_from) - psi(poss_pert_to_inds(poss_pert_to_range(1))-1,trans_from);
        lobs_prob_post = psi(poss_pert_to_inds(poss_pert_to_range(end)),trans_to) - psi(poss_pert_to_inds(poss_pert_to_range)-1,trans_to);
        lobs_prob = lobs_prob_pre + lobs_prob_post;
        for p = 1:length(poss_pert_to_range)
           prop_durs = cur_dur + delta_mat(:,poss_pert_to_range(p)); 
           poss_pert_from_range = find(prop_durs > 0 & prop_durs < size(P,2));
           log_duration_prob = p_d(prop_durs(poss_pert_from_range),trans_from);
           prob_llik = length_to(t-1,poss_pert_from_range) + log_duration_prob' + lobs_prob(p);
           if ~isempty(prob_llik)
               [cur_len,cur_loc] = max(prob_llik);
               length_to(t,p) = cur_len(1);
               cur_loc = cur_loc(1);
               coming_from(t,p) = poss_pert_from_range(cur_loc);
           else
               npdurs = size(delta_mat,1);
               length_to(t,p) = 1;
               coming_from(t,p) = round(npdurs/2);
           end
        end
    end
    
    [dummy,best_pert(num_trans)] = max(length_to(num_trans,:));
    for t=num_trans-1:-1:1,
        best_pert(t) = coming_from(t+1,best_pert(t+1));
    end
    new_state_trans = state_trans{i}(:,1) + pert_range(best_pert)';
    up_trans{i} = new_state_trans(state_trans{i}(:,2) == 2);
    down_trans{i} = new_state_trans(state_trans{i}(:,2) == 1);
    clear best_pert 
end

%% reconstruct optimized state sequence
for ns = 1:hmm.Nsegs
    new_hsmm_state_seq{ns} = ones(curT_new(ns),1);
    for i = 1:length(up_trans{ns})
        next_down = down_trans{ns}(find(down_trans{ns} > up_trans{ns}(i),1,'first'));
        if ~isempty(next_down)
            new_hsmm_state_seq{ns}(up_trans{ns}(i):next_down) = 2;
        else
            new_hsmm_state_seq{ns}(up_trans{ns}(i):end) = 2;
        end
    end
    if down_trans{ns}(1) < up_trans{ns}(1)
        new_hsmm_state_seq{ns}(1:down_trans{ns}(1)) = 2;
    end
end
