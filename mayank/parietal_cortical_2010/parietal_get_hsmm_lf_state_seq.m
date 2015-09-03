function [bb_hsmm_state_seq,hsmm_state_seq,hmm,Fs_bb,Fs] = parietal_get_hsmm_lf_state_seq(raw_data,raw_Fs,f_names)

addpath('F:\Code\smoothing\software')
addpath('F:\Code\FullBNT-1.0.4\KPMstats\')
addpath('F:\Code\FullBNT-1.0.4\netlab3.3')
addpath('F:\WC_Germany\new_stellate_analysis\')
addpath('F:\WC_Germany\hsmm_state_detection')

Fs_desired = 50; %(in Hz)
dsf = round(raw_Fs/Fs_desired);
dsf_bb = 8;
Fs_bb = raw_Fs/dsf_bb;

Fs = raw_Fs/dsf;
niqf = raw_Fs/2;
lf_lcf = 0.05; %low cut-off for the low-freq filter
lf_hcf = 2; %high cut-off for the low-freq filter
lf_hcf_d = 1; %high cut-off for the low-freq dist estimator
bb_lcf = 0.05;
bb_hcf = 20;

windowSize = 60; %window duration to estimate time varying state means (s)
windowSlide = 1; %window spacing to estimate time varying state means (s)
num_states = 2; %number of hidden states (must be 2 in current implementation)
state_mean_hcf = 0.01; %hcf for estimating time varying state means (0 for constant mean)

min_state_dur = 1;
max_state_dur = round(Fs*20);
dur_range = (1:max_state_dur)/Fs;
dur_hist_range = linspace(0,10,100);

%% initializations
%compute filter coefficients
[b_lf,a_lf] = butter(2,[lf_lcf/niqf lf_hcf/niqf]);
[b_bb,a_bb] = butter(2,[bb_lcf/niqf bb_hcf/niqf]);
[b_lfd,a_lfd] = butter(2,[lf_lcf/niqf lf_hcf_d/niqf]);
% [b_hf,a_hf] = butter(2,[hf_lcf/niqf hf_hcf/niqf]);

%compute low-freq observations
obs_lf = filtfilt(b_lf,a_lf,raw_data);
obs_lf = zscore(downsample(obs_lf,dsf));
obs_lf = obs_lf(:);

%compute low-freq observations
obs_bb = filtfilt(b_bb,a_bb,raw_data);
obs_bb = zscore(downsample(obs_bb,dsf_bb));
obs_bb = obs_bb(:);

obs_lfd = filtfilt(b_lfd,a_lfd,raw_data);
obs_lfd = zscore(downsample(obs_lfd,dsf));
obs_lfd = obs_lfd(:);

T = length(obs_lf);
t_axis = (1:T)/Fs;
t_axis_bb = (1:length(obs_bb))/Fs_bb;

%% initialize hsmm
emiss = obs_lfd;
[hmm] = initialize_hsmm(emiss,Fs,num_states,windowSize,windowSlide,state_mean_hcf);
emiss = obs_lf;

[hmm] = jmm_hmmtrain_hsmm(hmm,emiss);

[hmm_state_seq,lik_best,delta] = hsms_viterbi(hmm,emiss);

[smoothed_seq,orig_dur,reject_dur,sojourn_times] = thresh_state_smooth(hmm_state_seq,Fs,300,300);
[state_durations] = compute_state_durations(hmm_state_seq,Fs);
% n_hmm(1,:) = hist(state_durations{1},dur_hist_range);
% n_hmm(2,:) = hist(state_durations{2},dur_hist_range);
% n_hmm(1,:) = n_hmm(1,:)/sum(n_hmm(1,:));
% n_hmm(2,:) = n_hmm(2,:)/sum(n_hmm(2,:));

[state_durations] = compute_state_durations(smoothed_seq,Fs);
for i = 1:hmm.K
   if strcmp(hmm.state(i).dur_type,'inv_gauss')
       [params(1),params(2)] = inverse_gauss_mlfit(state_durations{i});
       hmm.state(i).dur_pars = params;
       hmm.P(i,:) = inverse_gaussian_pmf(dur_range,hmm.state(i).dur_pars(1),hmm.state(i).dur_pars(2));
   elseif strcmp(hmm.state(i).dur_type,'gamma')
       [params(1),params(2)] = gamma_mlfit(state_durations{i});
       hmm.state(i).dur_pars = params;
       hmm.P(i,:) = gamma_pmf(dur_range,hmm.state(i).dur_pars(1),hmm.state(i).dur_pars(2));
   else
       error('invalid duration distribution')
   end
end

%%
% hmm.state(1).alpha = alpha(1);
% hmm.state(1).beta = beta(1);
% hmm.state(2).alpha = alpha(2);
% hmm.state(2).beta = beta(2);

hmm.dur_range = dur_range;
hmm.max_state_dur = max_state_dur;
hmm.min_state_dur = min_state_dur;

% hmm.P(1,:) = gamma_pmf(dur_range,hmm.state(1).alpha,hmm.state(1).beta);
% hmm.P(2,:) = gamma_pmf(dur_range,hmm.state(2).alpha,hmm.state(2).beta);

% enforce any minimum state duration and renormalize
for i = 1:hmm.K
    if min_state_dur > 1
        hmm.P(i,1:min_state_dur-1) = zeros(1,min_state_dur-1);
        hmm.P(i,:) = hmm.P(i,:)/sum(hmm.P(i,:));
    end
end
    
numIt = 1;
[hmm2,gamma,hmm_window_post]=hsmm_train_real(hmm,emiss,numIt);

[hsmm_state_seq,max_lik] = hsms_viterbi_real(hmm2,emiss);

hmm = hmm2;
hmm.posteriors = hmm_window_post;
hmm.max_lik = max_lik;
hmm.gamma = gamma;

[bb_hsmm_state_seq] = pert_optimize_hsmm_transitions...
    (obs_bb,hmm,hsmm_state_seq,gamma,Fs,Fs_bb,'variable',0.15);


% [obs_dist,obsrange,win_t] = compute_slide_dist(obs_lf,windowSize,windowSlide,Fs);
% figure
% pcolor(win_t,obsrange,log(obs_dist'));shading flat
% caxis([-4 0])
% hold on
% plot(t_axis,hmm2.state(1).lf_meanfun,'w')
% plot(t_axis,hmm2.state(2).lf_meanfun,'w')
% t_names = ['G:\WC_Germany\parietal_cortical_2010\test_meanfuns\' f_names];
% print('-dpng',t_names), close


