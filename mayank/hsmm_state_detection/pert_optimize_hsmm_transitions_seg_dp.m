function [new_hsmm_state_seq] = pert_optimize_hsmm_transitions_seg_dp(new_obs,hmm,orig_state_seq,gamma,Fs_orig,Fs_new,pert_range)

%% currently works only for univariate nonstationary observations or
%% stationary multivariate observations

%these are the indices for the segments of data with UDS
UDS_samples = [];
for i = 1:hmm.Nsegs
    UDS_seg_inds{i} = hmm.UDS_segs(i,1):hmm.UDS_segs(i,2);
    curT(i) = length(UDS_seg_inds{i});
    UDS_samples = [UDS_samples hmm.UDS_segs(i,1):hmm.UDS_segs(i,2)];
    UDS_segind_new_start = round((hmm.UDS_segs(i,1)-1)*Fs_new/Fs_orig)+1;
    enddiff = length(new_obs)-UDS_segind_new_start+1;
    curT_new(i) = min(round(Fs_new/Fs_orig*curT(i)),enddiff);
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
new_obs = new_obs(:);

dsf = round(2016/Fs_orig);
dsf_new = round(2016/Fs_new);

[T_new,p] = size(new_obs);

A = hmm.A;

%compute state-dependent broadband amplitude
new_obs_d = downsample(new_obs,dsf/dsf_new);

if strcmp(meantype,'variable')
    meandiff = [];
    for i = 1:hmm.Nsegs
        t_axis{i} = (1:curT(i))/Fs_orig;
        t_axis_new{i} = (1:curT_new(i))/Fs_new;
        total_dur = curT(i)/hmm.Fs;
        numWins = floor((total_dur-hmm.windowSize)/hmm.windowSlide);
        win_t{i} = (0:numWins-1)*hmm.windowSlide+hmm.windowSize/2;
        win_t{i} = [0 win_t{i} max(t_axis_new{i})];
        temp_meanfun{i} = zeros(2,numWins);
        tot_prob = zeros(2,numWins);
        for l = 1:2
            for w = 1:numWins
                begT = round((w-1)*hmm.windowSlide*hmm.Fs)+1;
                endT = begT + round(hmm.windowSize*hmm.Fs);
                temp_meanfun{i}(l,w) = sum(gamma{i}(begT:endT,l).*new_obs_d(UDS_seg_inds{i}(begT:endT),:))/sum(gamma{i}(begT:endT,l));
                tot_prob(l,w) = sum(gamma{i}(begT:endT,l))/(hmm.windowSize*hmm.Fs);
            end
            temp_meanfun{i}(l,tot_prob(l,:) < 0.05) = nan;
        end
            meandiff = [meandiff temp_meanfun{i}(2,tot_prob(l,:) >= 0.05)-...
                temp_meanfun{i}(1,tot_prob(l,:) >= 0.05)];
    end
    avgmeandiff = nanmean(meandiff);
    for i = 1:hmm.Nsegs
        temp_meanfun{i}(1,isnan(temp_meanfun{i}(1,:))) = temp_meanfun{i}(2,isnan(temp_meanfun{i}(1,:))) - avgmeandiff;
        temp_meanfun{i}(2,isnan(temp_meanfun{i}(2,:))) = temp_meanfun{i}(1,isnan(temp_meanfun{i}(2,:))) + avgmeandiff;
        temp_meanfun{i} = [temp_meanfun{i}(:,1) temp_meanfun{i} temp_meanfun{i}(:,end)]';
        meanfun_est{i} = interp1(win_t{i},temp_meanfun{i},t_axis_new{i});
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
           mdiff{i} = new_obs(UDS_seg_inds_new{i},:) - meanfun_est{i}(:,l);
           o_sum = o_sum + (gamma_inter{i}(:,l).*mdiff{i})'*mdiff{i};
           g_sum = g_sum + sum(gamma_inter{i}(:,l));
        end
        var_est(l) = o_sum/g_sum;
        for i = 1:hmm.Nsegs
            state_llik{i}(l,:) = -mdiff{i}.^2/(2*var_est(l));
            state_llik{i}(l,:) = state_llik{i}(l,:) - 1/2*log(2*pi*var_est(l));
        end
    end
elseif strcmp(meantype,'fixed')
    for i = 1:hmm.Nsegs
        t_axis{i}(1) = 0;
        t_axis{i}(end) = t_axis_new{i}(end);
        gamma_inter{i} = interp1(t_axis{i},gamma{i},t_axis_new{i});
    end
    for l = 1:hmm.K
        o_sum = 0;
        g_sum = 0;
        for i = 1:hmm.Nsegs
            o_sum = o_sum + repmat(gamma_inter{i}(:,l),p,1).*new_obs(UDS_seg_inds_new{i},:);
            g_sum = g_sum + sum(gamma_inter{i}(:,l));
        end
        mean_est = o_sum/g_sum;
        o_sum = 0;
        for i = 1:hmm.Nsegs
           mdiff{i} = new_obs(UDS_seg_inds_new{i},:) - repmat(mean_est',curT_new(i),1);          
           o_sum = o_sum + (repmat(gamma_inter{i}(:,l),1,p).*mdiff{i})'*mdiff{i};
        end
        var_est = o_sum/g_sum;
        [T,err] = cholcov(var_est,0);
        if err ~= 0
            error('Invalid covariance matrix')
        end
        for i = 1:hmm.Nsegs
            X0 = mdiff{i}/T;
            quadform = sum(X0.^2,2);
            state_llik{i}(l,:) = -1/2*log((2*pi)^p*det(var_est)) - 1/2*quadform;
        end
    end
else
    error('invalid mean type')
end


for i = 1:hmm.Nsegs
    obs_eps = -200;
    state_llik{i}(state_llik{i} < obs_eps) = obs_eps;
    %find instances where the emissions model fails and correct the
    %likelihood
    false_ups = find(new_obs(UDS_seg_inds_new{i},:)' < meanfun_est{i}(:,1)' & ...
        state_llik{i}(1,:) < state_llik{i}(2,:));
    state_llik{i}(1,false_ups) = state_llik{i}(2,false_ups)+1;
    false_downs = find(new_obs(UDS_seg_inds_new{i},:)' > meanfun_est{i}(:,2)' & ...
        state_llik{i}(1,:) > state_llik{i}(2,:));
     state_llik{i}(2,false_downs) = state_llik{i}(1,false_downs)+1;
end

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
    else
        P(k,:) = inverse_gaussian_pmf(poss_durs/Fs_new,hmm.state(k).dur_pars(1),hmm.state(k).dur_pars(2));
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
           [cur_len,cur_loc] = max(prob_llik);
           length_to(t,p) = cur_len(1);
           cur_loc = cur_loc(1);
           coming_from(t,p) = poss_pert_from_range(cur_loc);
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
