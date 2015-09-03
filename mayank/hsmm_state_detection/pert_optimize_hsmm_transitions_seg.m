function [new_hsmm_state_seq,no_peak_ups,no_peak_downs] = pert_optimize_hsmm_transitions_seg(new_obs,hmm,orig_state_seq,gamma,Fs_orig,Fs_new,up_max_pert,down_max_pert)

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


%% for each detected up transition find most likely transition time from state posteriors
if length(up_max_pert) > 1
    up_back_pert = round(Fs_new*up_max_pert(1));
    up_forward_pert = round(Fs_new*up_max_pert(2));
else
    up_back_pert = round(Fs_new*up_max_pert);
    up_forward_pert = round(Fs_new*up_max_pert);
end
if length(down_max_pert) > 1
    down_back_pert = round(Fs_new*down_max_pert(1));
    down_forward_pert = round(Fs_new*down_max_pert(2));
else
    down_back_pert = round(Fs_new*down_max_pert);
    down_forward_pert = round(Fs_new*down_max_pert);
end
for i = 1:hmm.Nsegs
    up_trans{i} = find(orig_state_seq{i}(1:end-1) == 1 & orig_state_seq{i}(2:end) == 2);
    down_trans{i} = find(orig_state_seq{i}(1:end-1) == 2 & orig_state_seq{i}(2:end) == 1);
    up_trans{i} = round(up_trans{i}*Fs_new/Fs_orig);
    down_trans{i} = round(down_trans{i}*Fs_new/Fs_orig);
    up_trans{i}(up_trans{i} <= max(up_back_pert,down_back_pert) | up_trans{i} + max(up_forward_pert,down_forward_pert) >= curT_new(i)) = [];
    down_trans{i}(down_trans{i} <= max(up_back_pert,down_back_pert) | down_trans{i} + max(up_forward_pert,down_forward_pert) >= curT_new(i)) = [];
end


%% optimize up transitions
% if trans_flag(1)
for ns = 1:hmm.Nsegs
    no_peak_ups{ns} = [];
    for i = 1:length(up_trans{ns})
        cur_up = up_trans{ns}(i);
        prev_down = down_trans{ns}(find(down_trans{ns} < cur_up,1,'last'));
        if ~isempty(prev_down)
            prev_dur = (cur_up - prev_down)/Fs_new;
        else
            prev_dur = cur_up/Fs_new;
        end
        next_down = down_trans{ns}(find(down_trans{ns} > cur_up,1,'first'));
        if ~isempty(next_down)
            next_dur = (next_down - cur_up)/Fs_new;
        else
            next_dur = (curT_new(ns)-cur_up)/Fs_new;
        end
        if strcmp(hmm.state(1).dur_type,'gamma')
            prev_state_probs = log(gamma_pmf(prev_dur + (-up_back_pert:up_forward_pert)/Fs_new,hmm.state(1).dur_pars(1),hmm.state(1).dur_pars(2)));
        elseif strcmp(hmm.state(1).dur_type,'inv_gauss')
            prev_state_probs = log(inverse_gaussian_pmf(prev_dur + (-up_back_pert:up_forward_pert)/Fs_new,hmm.state(1).dur_pars(1),hmm.state(1).dur_pars(2)));
        else
            prev_state_probs = log(A(1,1).^(round(prev_dur*Fs_new) + (-up_back_pert:up_forward_pert))*A(1,2));
        end
        if strcmp(hmm.state(2).dur_type,'gamma')
            next_state_probs = log(gamma_pmf(next_dur - (-up_back_pert:up_forward_pert)/Fs_new,hmm.state(2).dur_pars(1),hmm.state(2).dur_pars(2)));
        elseif strcmp(hmm.state(2).dur_type,'inv_gauss')
            next_state_probs = log(inverse_gaussian_pmf(next_dur - (-up_back_pert:up_forward_pert)/Fs_new,hmm.state(2).dur_pars(1),hmm.state(2).dur_pars(2)));
        else
            next_state_probs = log(A(2,2).^(round(next_dur*Fs_new) - (-up_back_pert:up_forward_pert))*A(2,1));
        end
        poss_ups = cur_up - up_back_pert:cur_up + up_forward_pert;
        cur_state_llik = state_llik{ns}(:,poss_ups);
        down_lliks = cumsum(cur_state_llik(1,:));
        up_lliks = fliplr(cumsum(fliplr(cur_state_llik(2,:))));
        [dummy,peakloc] = max(down_lliks+up_lliks+prev_state_probs+next_state_probs);
        [dummy,peaklocs] = findpeaks(down_lliks+up_lliks);
        if ~isempty(peaklocs)
            [dummy,bestpeak] = min(abs(cur_up - peaklocs));
            peakloc = peaklocs(bestpeak);
            up_trans{ns}(i) = peakloc + cur_up - up_back_pert;
        else
            no_peak_ups{ns} = [no_peak_ups{ns} i];
        end
    end
    % end
end

%% optimize down transitions
% if trans_flag(2)
for ns = 1:hmm.Nsegs
    no_peak_downs{ns} = [];
    for i = 1:length(down_trans{ns})
        cur_down = down_trans{ns}(i);
        prev_up = up_trans{ns}(find(up_trans{ns} < cur_down,1,'last'));
        if ~isempty(prev_up)
            prev_dur = (cur_down - prev_up)/Fs_new;
        else
            prev_dur = cur_down/Fs_new;
        end
        next_up = up_trans{ns}(find(up_trans{ns} > cur_down,1,'first'));
        if ~isempty(next_up)
            next_dur = (next_up - cur_down)/Fs_new;
        else
            next_dur = (curT_new(ns)-cur_down)/Fs_new;
        end
        if strcmp(hmm.state(2).dur_type,'gamma')
            prev_state_probs = log(gamma_pmf(prev_dur + (-down_back_pert:down_forward_pert)/Fs_new,hmm.state(2).dur_pars(1),hmm.state(2).dur_pars(2)));
        elseif strcmp(hmm.state(2).dur_type,'inv_gauss')
            prev_state_probs = log(inverse_gaussian_pmf(prev_dur + (-down_back_pert:down_forward_pert)/Fs_new,hmm.state(2).dur_pars(1),hmm.state(2).dur_pars(2)));
        else
            prev_state_probs = log(A(2,2).^(round(prev_dur*Fs_new) + (-down_back_pert:down_forward_pert))*A(2,1));
        end
        if strcmp(hmm.state(1).dur_type,'gamma')
            next_state_probs = log(gamma_pmf(next_dur - (-down_back_pert:down_forward_pert)/Fs_new,hmm.state(1).dur_pars(1),hmm.state(1).dur_pars(2)));
        elseif strcmp(hmm.state(1).dur_type,'inv_gauss')
            next_state_probs = log(inverse_gaussian_pmf(next_dur - (-down_back_pert:down_forward_pert)/Fs_new,hmm.state(1).dur_pars(1),hmm.state(1).dur_pars(2)));
        else
            next_state_probs = log(A(1,1).^(round(next_dur*Fs_new) - (-down_back_pert:down_forward_pert))*A(1,2));
        end
        poss_downs = cur_down - down_back_pert:cur_down + down_forward_pert;
        cur_state_llik = state_llik{ns}(:,poss_downs);
        up_lliks = cumsum(cur_state_llik(2,:));
        down_lliks = fliplr(cumsum(fliplr(cur_state_llik(1,:))));
        [dummy,peaklocs] = findpeaks(down_lliks+up_lliks);
        if ~isempty(peaklocs)
            [dummy,bestpeak] = min(abs(cur_down - peaklocs));
            peakloc = peaklocs(bestpeak);
            down_trans{ns}(i) = peakloc + cur_down - down_back_pert;
        else
            no_peak_downs{ns} = [no_peak_downs{ns} i];
        end
    end
% end
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
