function [bb_hsmm_state_seq,hsmm_state_seq,bb_hmm_state_seq,hsmm,hmm,Fs_bb,Fs] = ...
    parietal_get_hsmm_uds_state_seq_fixmean_new(raw_data,raw_Fs,desynch_times,f_names)

tic;
drive_letter = 'G';
addpath(strcat(drive_letter,':\WC_Germany\hsmm_state_detection'))

Fs_desired = 50; %(in Hz)
dsf = round(raw_Fs/Fs_desired);

Fs = raw_Fs/dsf;
niqf = raw_Fs/2;
lf_lcf = 0.05; %low cut-off for the low-freq filter
lf_hcf = 2; %high cut-off for the low-freq filter
num_states = 2; %number of hidden states (must be 2 in current implementation)

min_state_dur = 1;
max_state_dur = round(Fs*30);
min_seg_dur = 60;
dur_range = (1:max_state_dur)/Fs;

%% initializations
[obs_lf,t_axis,Fs] = get_lf_features(raw_data,raw_Fs,Fs_desired,[lf_lcf lf_hcf]);

%% compute UDS segments
UDS_segs = hsmm_uds_get_uds_segments(desynch_times,Fs,length(obs_lf),min_seg_dur);

%% initialize hsmm
emiss = obs_lf;
params.meantype = 'fixed';
params.UDS_segs = UDS_segs;
params.movingwin = [50 10];
[hmm] = hsmm_uds_initialize(emiss,Fs,params);

%%
hmm.min_mllik = -8; %set minimum value of the maximal log likelihood for either state before switching to robust estimator.  If you don't want to use this feature, set this to -Inf
hmm.windowSize = 50;
hmm = hsmm_uds_train_hmm(hmm,emiss);
hsmm = hmm;
%%
[hmm_state_seq,llik_best] = hsmm_uds_viterbi_hmm(hmm,emiss);

hmm.run_time = toc;

smoothed_seq = thresh_state_smooth_seg(hmm_state_seq,Fs,100,100); %smooth out excessively short state durations to improve model fits
[state_durations] = compute_state_durations_seg(smoothed_seq,Fs); %compute the set of up and down state durations
for i = 1:hmm.K
    %estimate the empirical state duration pmf
    emp_pmf = hist(state_durations{i},dur_range);
    emp_pmf = emp_pmf/sum(emp_pmf);
    [mu,lambda] = inverse_gauss_mlfit(state_durations{i});%compute the ML parameters of the inverse gaussian dist
    ig_pmf = inverse_gaussian_pmf(dur_range,mu,lambda);%estimate the IG PMF
    [alpha,beta] = gamma_mlfit(state_durations{i}); %compute the ML parameters of the gamma
    gam_pmf = gamma_pmf(dur_range,alpha,beta); %estimate the gamma PMF
    ks_ig = max(abs(cumsum(emp_pmf) - cumsum(ig_pmf))); %KS stat for the best-fitting IG dist
    ks_gam = max(abs(cumsum(emp_pmf) - cumsum(gam_pmf))); %KS stat for the best-fitting GAM dist
    hsmm.state(i).dur_type = 'inv_gauss'; %use inverse gaussian model in all cases
    hsmm.state(i).dur_pars = [mu lambda];
    hsmm.P(i,:) = ig_pmf;
end

%%
%initialize duration model parameters
hsmm.dur_range = dur_range;
hsmm.max_state_dur = max_state_dur;
hsmm.min_state_dur = min_state_dur;
hmm.dur_range = dur_range;
hmm.max_state_dur = max_state_dur;
hmm.min_state_dur = min_state_dur;

% enforce any minimum state duration and renormalize
for i = 1:hmm.K
    if min_state_dur > 1
        hsmm.P(i,1:min_state_dur-1) = zeros(1,min_state_dur-1);
        hsmm.P(i,:) = hsmm.P(i,:)/sum(hsmm.P(i,:));
    end
end
    
% [hmm,gamma,hmm_window_post,bad_model]=hsmm_uds_train_seg(hmm,emiss);
[hsmm,gamma,hmm_window_post]=hsmm_uds_train_hsmm(hsmm,emiss);

% [hsmm_state_seq,max_lik] = hsmm_uds_viterbi_seg(hmm,emiss);
[hsmm_state_seq,max_lik] = hsmm_uds_viterbi_hsmm(hsmm,emiss);

hsmm.posteriors = hmm_window_post;
hsmm.max_lik = max_lik;

%% for optimizing transition times
dsf_bb = 8;
Fs_bb = raw_Fs/dsf_bb;
bb_lcf = 0.05;
bb_hcf = 20;
[obs_bb,t_axis,Fs_bb] = get_lf_features(raw_data,raw_Fs,Fs_bb,[bb_lcf bb_hcf]);

pert_range_bb = [-0.15 0.15];
% [bb_hsmm_state_seq] = pert_optimize_hsmm_transitions_seg_dp_robust... 
%     (obs_bb,hsmm,hsmm_state_seq,gamma,Fs,Fs_bb,pert_range_bb);
[bb_hsmm_state_seq] = hsmm_uds_pert_optimize_transitions_v2(obs_bb,hsmm,hsmm_state_seq,gamma,Fs,Fs_bb,pert_range_bb);

hsmm.run_time = toc;

temp_p = zeros(size(hsmm.P));
temp_p(1,:) = geometric_pmf(1:length(dur_range),hmm.A(2,1));
temp_p(2,:) = geometric_pmf(1:length(dur_range),hmm.A(1,2));
hmm.P = temp_p;
for i = 1:hmm.K
    hmm.state(i).dur_type = 'hmm';
end

% [bb_hmm_state_seq] = pert_optimize_hsmm_transitions_seg_dp_robust... 
%     (obs_bb,hmm_temp,hmm_state_seq,gamma_hmm,Fs,Fs_bb,pert_range_bb);
[bb_hmm_state_seq] = hsmm_uds_pert_optimize_transitions_v2(obs_bb,hmm,hmm_state_seq,hmm.gamma,Fs,Fs_bb,pert_range_bb);


% [bb_hmm_state_seq] = pert_optimize_hsmm_transitions_seg...
%     (obs_bb,hmm,hmm_state_seq,gamma_hmm,Fs,Fs_bb,up_max_pertbb,down_max_pertbb);
% [hf_hsmm_state_seq] = pert_optimize_hsmm_transitions_seg_dp_robust...
%     (obs_hf,hmm,hsmm_state_seq,gamma,Fs,Fs_bb,up_max_perthf,down_max_perthf);


% [obs_dist,obsrange,win_t] = compute_slide_dist(obs_lf,meanparams.windowSize,meanparams.windowSlide,Fs);
% figure
% pcolor(win_t,obsrange,log(obs_dist'));shading flat
% caxis([-4 0])
% hold on
% for i = 1:hmm.Nsegs
%     t_axis = (hmm.UDS_segs(i,1):hmm.UDS_segs(i,2))/hmm.Fs;
%     plot(t_axis,hmm.state(1).meanfun{i},'w')
%     plot(t_axis,hmm.state(2).meanfun{i},'w')
% end
% t_names = ['F:\WC_Germany\parietal_cortical_2010\state_meanfuns\' f_names];
% print('-dpng',t_names), close


