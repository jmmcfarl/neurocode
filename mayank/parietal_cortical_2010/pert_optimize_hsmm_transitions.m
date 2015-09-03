function [new_hsmm_state_seq,no_peak_ups,no_peak_downs] = pert_optimize_hsmm_transitions(new_obs,hmm,orig_state_seq,gamma,Fs_orig,Fs_new,max_pert)

%% currently works only for univariate nonstationary observations or
%% stationary multivariate observations

meantype = hmm.meantype;
if size(gamma,1) < size(gamma,2)
    gamma= gamma';
end
new_obs = new_obs(:);

no_peak_ups = [];
no_peak_downs = [];

dsf = round(2016/Fs_orig);
dsf_new = round(2016/Fs_new);

[T_new,p] = size(new_obs);

A = hmm.A;

t_axis = (1:hmm.T)/Fs_orig;
t_axis_new = (1:T_new)/Fs_new;

%compute state-dependent broadband amplitude
new_obs_d = downsample(new_obs,dsf/dsf_new);

if strcmp(meantype,'variable')
    total_dur = hmm.T/hmm.Fs;
    numWins = floor((total_dur-hmm.windowSize)/hmm.windowSlide);
    win_t = (0:numWins-1)*hmm.windowSlide+hmm.windowSize/2;
    win_t = [0 win_t max(t_axis_new)];
    temp_meanfun = zeros(2,numWins);
    tot_prob = zeros(2,numWins);
    for l = 1:2
        for w = 1:numWins
            begT = round((w-1)*hmm.windowSlide*hmm.Fs)+1;
            endT = begT + round(hmm.windowSize*hmm.Fs);
            temp_meanfun(l,w) = sum(gamma(begT:endT,l).*new_obs_d(begT:endT,:))/sum(gamma(begT:endT,l));
            tot_prob(l,w) = sum(gamma(begT:endT,l))/(hmm.windowSize*hmm.Fs);
        end
        temp_meanfun(l,tot_prob(l,:) < 0.05) = nan;
    end
    meandiff = temp_meanfun(2,:)-temp_meanfun(1,:);
    meandiff = nanmean(meandiff);
    temp_meanfun(1,isnan(temp_meanfun(1,:))) = temp_meanfun(2,isnan(temp_meanfun(1,:))) - meandiff;
    temp_meanfun(2,isnan(temp_meanfun(2,:))) = temp_meanfun(1,isnan(temp_meanfun(2,:))) + meandiff;
    temp_meanfun = [temp_meanfun(:,1) temp_meanfun temp_meanfun(:,end)]';
    meanfun_est = interp1(win_t,temp_meanfun,t_axis_new);
    t_axis(1) = 0;
    t_axis(end) = t_axis_new(end);
    gamma_inter = interp1(t_axis,gamma,t_axis_new);
    for i = 1:hmm.K
        mdiff = new_obs-meanfun_est(:,i);
        var_est(i) = ((gamma_inter(:,i).*mdiff)'*mdiff)/sum(gamma_inter(:,i));
        state_llik(i,:) = -mdiff.^2/(2*var_est(i));
        state_llik(i,:) = state_llik(i,:) - 1/2*log(2*pi*var_est(i));
    end
elseif strcmp(meantype,'fixed')
    t_axis(1) = 0;
    t_axis(end) = t_axis_new(end);
    gamma_inter = interp1(t_axis,gamma,t_axis_new);
    for i = 1:hmm.K
        mean_est = sum(repmat(gamma_inter(:,i),p,1).*new_obs)/sum(gamma_inter(:,i));
        mdiff = new_obs - repmat(mean_est',T_new,1);
        var_est = ((repmat(gamma_inter(:,i),1,p).*mdiff)'*mdiff)/sum(gamma_inter(:,i));
        [T,err] = cholcov(var_est,0);
        if err ~= 0
            error('Invalid covariance matrix')
        end
        X0 = mdiff/T;
        quadform = sum(X0.^2,2);
        state_llik(i,:) = -1/2*log((2*pi)^p*det(var_est)) - 1/2*quadform;
    end
else
    error('invalid mean type')
end


%% for each detected up transition find most likely transition time from state posteriors
back_pert = round(Fs_new*max_pert);
forward_pert = round(Fs_new*max_pert);
up_trans = find(orig_state_seq(1:end-1) == 1 & orig_state_seq(2:end) == 2);
down_trans = find(orig_state_seq(1:end-1) == 2 & orig_state_seq(2:end) == 1);
up_trans = round(up_trans*Fs_new/Fs_orig);
down_trans = round(down_trans*Fs_new/Fs_orig);
up_trans(up_trans <= back_pert | up_trans + forward_pert >= length(new_obs)) = [];
down_trans(down_trans <= back_pert | down_trans + forward_pert >= length(new_obs)) = [];

%% optimize up transitions
% if trans_flag(1)
    for i = 1:length(up_trans)
        cur_up = up_trans(i);
        prev_down = down_trans(find(down_trans < cur_up,1,'last'));
        if ~isempty(prev_down)
            prev_dur = (cur_up - prev_down)/Fs_new;
        else
            prev_dur = cur_up/Fs_new;
        end
        next_down = down_trans(find(down_trans > cur_up,1,'first'));
        if ~isempty(next_down)
            next_dur = (next_down - cur_up)/Fs_new;
        else
            next_dur = (T_new-cur_up)/Fs_new;
        end
        if strcmp(hmm.state(1).dur_type,'gamma')
            prev_state_probs = log(gamma_pmf(prev_dur + (-back_pert:forward_pert)/Fs_new,hmm.state(1).dur_pars(1),hmm.state(1).dur_pars(2)));
        elseif strcmp(hmm.state(1).dur_type,'inv_gauss')
            prev_state_probs = log(inverse_gaussian_pmf(prev_dur + (-back_pert:forward_pert)/Fs_new,hmm.state(1).dur_pars(1),hmm.state(1).dur_pars(2)));
        else
            prev_state_probs = log(A(1,1).^(round(prev_dur*Fs_new) + (-back_pert:forward_pert))*A(1,2));
        end
        if strcmp(hmm.state(2).dur_type,'gamma')
            next_state_probs = log(gamma_pmf(next_dur - (-back_pert:forward_pert)/Fs_new,hmm.state(2).dur_pars(1),hmm.state(2).dur_pars(2)));
        elseif strcmp(hmm.state(2).dur_type,'inv_gauss')
            next_state_probs = log(inverse_gaussian_pmf(next_dur - (-back_pert:forward_pert)/Fs_new,hmm.state(2).dur_pars(1),hmm.state(2).dur_pars(2)));
        else
            next_state_probs = log(A(2,2).^(round(next_dur*Fs_new) - (-back_pert:forward_pert))*A(2,1));
        end
        poss_ups = cur_up - back_pert:cur_up + forward_pert;
        cur_state_llik = state_llik(:,poss_ups);
        down_lliks = cumsum(cur_state_llik(1,:));
        up_lliks = fliplr(cumsum(fliplr(cur_state_llik(2,:))));
        [dummy,peakloc] = max(down_lliks+up_lliks+prev_state_probs+next_state_probs);
        [dummy,peaklocs] = findpeaks(down_lliks+up_lliks);
        if ~isempty(peaklocs)
            [dummy,bestpeak] = min(abs(cur_up - peaklocs));
            peakloc = peaklocs(bestpeak);
            up_trans(i) = peakloc + cur_up - back_pert;
        else
            no_peak_ups = [no_peak_ups i];
        end
    end
% end

%% optimize down transitions
% if trans_flag(2)
    for i = 1:length(down_trans)
        cur_down = down_trans(i);
        prev_up = up_trans(find(up_trans < cur_down,1,'last'));
        if ~isempty(prev_up)
            prev_dur = (cur_down - prev_up)/Fs_new;
        else
            prev_dur = cur_down/Fs_new;
        end
        next_up = up_trans(find(up_trans > cur_down,1,'first'));
        if ~isempty(next_up)
            next_dur = (next_up - cur_down)/Fs_new;
        else
            next_dur = (length(new_obs)-cur_down)/Fs_new;
        end
        if strcmp(hmm.state(2).dur_type,'gamma')
            prev_state_probs = log(gamma_pmf(prev_dur + (-back_pert:forward_pert)/Fs_new,hmm.state(2).dur_pars(1),hmm.state(2).dur_pars(2)));
        elseif strcmp(hmm.state(2).dur_type,'inv_gauss')
            prev_state_probs = log(inverse_gaussian_pmf(prev_dur + (-back_pert:forward_pert)/Fs_new,hmm.state(2).dur_pars(1),hmm.state(2).dur_pars(2)));
        else
            prev_state_probs = log(A(2,2).^(round(prev_dur*Fs_new) + (-back_pert:forward_pert))*A(2,1));
        end
        if strcmp(hmm.state(1).dur_type,'gamma')
            next_state_probs = log(gamma_pmf(next_dur - (-back_pert:forward_pert)/Fs_new,hmm.state(1).dur_pars(1),hmm.state(1).dur_pars(2)));
        elseif strcmp(hmm.state(1).dur_type,'inv_gauss')
            next_state_probs = log(inverse_gaussian_pmf(next_dur - (-back_pert:forward_pert)/Fs_new,hmm.state(1).dur_pars(1),hmm.state(1).dur_pars(2)));
        else
            next_state_probs = log(A(1,1).^(round(next_dur*Fs_new) - (-back_pert:forward_pert))*A(1,2));
        end
        poss_downs = cur_down - back_pert:cur_down + forward_pert;
        cur_state_llik = state_llik(:,poss_downs);
        up_lliks = cumsum(cur_state_llik(2,:));
        down_lliks = fliplr(cumsum(fliplr(cur_state_llik(1,:))));
        [dummy,peaklocs] = findpeaks(down_lliks+up_lliks);
        if ~isempty(peaklocs)
            [dummy,bestpeak] = min(abs(cur_down - peaklocs));
            peakloc = peaklocs(bestpeak);
            down_trans(i) = peakloc + cur_down - back_pert;
        else
            no_peak_downs = [no_peak_downs i];
        end
    end
% end

%% reconstruct optimized state sequence
new_hsmm_state_seq = ones(T_new,1);
for i = 1:length(up_trans)
    next_down = down_trans(find(down_trans > up_trans(i),1,'first'));
    if ~isempty(next_down)
        new_hsmm_state_seq(up_trans(i):next_down) = 2;
    else
        new_hsmm_state_seq(up_trans(i):end) = 2;
    end
end
if down_trans(1) < up_trans(1)
    new_hsmm_state_seq(1:down_trans(1)) = 2;
end

