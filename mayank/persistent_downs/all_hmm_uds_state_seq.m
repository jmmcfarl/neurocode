function [bb_hmm_state_seq,hmm,Fs_bb,Fs] = all_hmm_uds_state_seq(raw_data,raw_Fs,desynch_times,f_names)

drive_letter = 'C';
addpath(strcat(drive_letter,':\WC_Germany\hsmm_state_detection'))

Fs_desired = 50; %(in Hz)
dsf = round(raw_Fs/Fs_desired);

Fs = raw_Fs/dsf;
niqf = raw_Fs/2;
lf_lcf = 0.05; %low cut-off for the low-freq filter (default 0.05)
lf_hcf = 2; %high cut-off for the low-freq filter
num_states = 2; %number of hidden states (must be 2 in current implementation)

meanparams.windowSize = 50; %window duration to estimate time varying state means (s)
meanparams.windowSlide = 10; %window spacing to estimate time varying state means (s)
meanparams.meantype = 'variable';

min_state_dur = 1;
max_state_dur = round(Fs*60);
min_seg_dur = 60;
dur_range = (1:max_state_dur)/Fs;

%% initializations
[obs_lf,t_axis,Fs] = get_lf_features(raw_data,raw_Fs,Fs_desired,[lf_lcf lf_hcf]);

%% compute UDS segments
UDS_segs = hsmm_uds_get_uds_segments(desynch_times,Fs,length(obs_lf),min_seg_dur);

%% initialize hsmm
emiss = obs_lf;
params.meantype = 'variable';
params.UDS_segs = UDS_segs;
params.movingwin = [50 10];
[hmm] = hsmm_uds_initialize(emiss,Fs,params);

%%
hmm.min_mllik = -8; %set minimum value of the maximal log likelihood for either state before switching to robust estimator.  If you don't want to use this feature, set this to -Inf
hmm.windowSize = 50; %default 50
hmm = hsmm_uds_train_hmm(hmm,emiss);

%%
[hmm_state_seq,llik_best] = hsmm_uds_viterbi_hmm(hmm,emiss);

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
end

%% for optimizing transition times
dsf_bb = 8;
Fs_bb = raw_Fs/dsf_bb;
bb_lcf = 0.05;
bb_hcf = 20;
[obs_bb,t_axis,Fs_bb] = get_lf_features(raw_data,raw_Fs,Fs_bb,[bb_lcf bb_hcf]);

pert_range_bb = [-0.15 0.15];
temp_p = zeros(2,length(dur_range));
temp_p(1,:) = geometric_pmf(1:length(dur_range),hmm.A(2,1));
temp_p(2,:) = geometric_pmf(1:length(dur_range),hmm.A(1,2));
hmm.P = temp_p;
for i = 1:hmm.K
    hmm.state(i).dur_type = 'hmm';
end
hmm.dur_range = dur_range;
hmm.max_state_dur = max_state_dur;
hmm.min_state_dur = min_state_dur;
[bb_hmm_state_seq] = hsmm_uds_pert_optimize_transitions_v2(obs_bb,hmm,hmm_state_seq,hmm.gamma,Fs,Fs_bb,pert_range_bb);
[obs_dist,obsrange,win_t] = compute_slide_dist(obs_lf,meanparams.windowSize,meanparams.windowSlide,Fs);

figure('visible','off')
pcolor(win_t,obsrange,log(obs_dist'));shading flat
caxis([-4 0])
hold on
for i = 1:hmm.Nsegs
    t_axis = (hmm.UDS_segs(i,1):hmm.UDS_segs(i,2))/hmm.Fs;
    plot(t_axis,hmm.state(1).meanfun{i},'w')
    plot(t_axis,hmm.state(2).meanfun{i},'w')
end
t_names = ['C:\WC_Germany\overall_EC\new_state_meanfuns\' f_names];
% print('-dpng',t_names)
close


