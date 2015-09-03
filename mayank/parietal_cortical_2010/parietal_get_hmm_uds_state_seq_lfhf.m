function [bb_hmm_state_seq,hmm,Fs_bb,Fs] = ...
    parietal_get_hmm_uds_state_seq_lfhf(raw_data,raw_Fs,desynch_times)

drive_letter = 'G';
addpath(strcat(drive_letter,':\Code\FullBNT-1.0.4\KPMstats\'))
addpath(strcat(drive_letter,':\Code\FullBNT-1.0.4\netlab3.3'))
addpath(strcat(drive_letter,':\WC_Germany\hsmm_state_detection'))

Fs_desired = 50; %(in Hz)
dsf = round(raw_Fs/Fs_desired);

Fs = raw_Fs/dsf;
niqf = raw_Fs/2;
lf_lcf = 0.05; %low cut-off for the low-freq filter
lf_hcf = 2; %high cut-off for the low-freq filter
lf_hcf_d = 1; %high cut-off for the low-freq dist estimator
hf_lcf = 20; %low cut-off for the low-freq filter
hf_hcf = 80; %high cut-off for the low-freq filter
hf_smooth = 0.15;
num_states = 2; %number of hidden states (must be 2 in current implementation)

meanparams.windowSize = 50; %window duration to estimate time varying state means (s)
meanparams.windowSlide = 1; %window spacing to estimate time varying state means (s)
% meanparams.state_mean_hcf = 0.01; %hcf for estimating time varying state means (0 for constant mean)
meanparams.meantype = 'variable';

min_state_dur = 1;
max_state_dur = round(Fs*30);
min_seg_dur = 60;
dur_range = (1:max_state_dur)/Fs;
% dur_hist_range = linspace(0,10,100);

%% initializations
[obs_lf,t_axis,Fs] = get_lf_features(raw_data,raw_Fs,Fs_desired,[lf_lcf lf_hcf]);
[obs_hf,t_axis,Fs] = get_hf_features(raw_data,raw_Fs,Fs_desired,[hf_lcf hf_hcf],hf_smooth);
% [obs_lfd,t_axis,Fs] = get_lf_features(raw_data,raw_Fs,Fs_desired,[lf_lcf lf_hcf_d]);

%% compute UDS segments
N = length(obs_lf);
desynch_indices = round(desynch_times*Fs);
temp_inds = zeros(N,1);
for i = 1:size(desynch_indices,1)
    temp_inds(desynch_indices(i,1):desynch_indices(i,2)) = 1;
end
UDS_start = find(temp_inds(1:end-1)==1 & temp_inds(2:end)==0);
UDS_stop = find(temp_inds(1:end-1)==0 & temp_inds(2:end)==1);
if temp_inds(1)==0
    UDS_start = [1; UDS_start];
end
if temp_inds(end)==0
    UDS_stop = [UDS_stop; N];
end
UDS_segs = [UDS_start UDS_stop];

UDS_seg_durs = (UDS_segs(:,2)-UDS_segs(:,1))/Fs;
too_short = find(UDS_seg_durs < min_seg_dur);
UDS_segs(too_short,:) = [];

%% initialize hsmm
% [hmm_lf] = initialize_uds_hsmm_seg(obs_lfd,Fs,meanparams,UDS_segs,drive_letter);
% [hmm_hf] = initialize_uds_hsmm_seg(obs_hf,Fs,meanparams,UDS_segs,drive_letter);
% emiss = [obs_lfd obs_hf];
emiss = [obs_lf obs_hf];
hmm = initialize_uds_hsmm_seg_dualvar(emiss,Fs,meanparams,UDS_segs,drive_letter);
% emiss(:,1) = obs_lf;
hmm.min_mllik = -3; %set minimum value of the maximal log likelihood for either state before switching to robust estimator

% hmm = hmm_lf;
% hmm.p = 2;
% hmm.state(1).meanfun2 = hmm_hf.state(1).meanfun;
% hmm.state(2).meanfun2 = hmm_hf.state(2).meanfun;
% ov_cov = cov(obs_lf,obs_hf);
% ov_cov = ov_cov(1,2);
% hmm.state(1).var = [hmm.state(1).var ov_cov*0.5;
%     ov_cov*0.5 hmm_hf.state(1).var];
% hmm.state(2).var = [hmm.state(2).var ov_cov*hmm.state(2).var;
%     ov_cov*hmm.state(2).var hmm_hf.state(2).var];

%%
[hmm, gamma_hmm] = hmm_uds_train_seg_duallfp(hmm,emiss);

[hmm_state_seq,llik_best] = hmm_uds_viterbi_seg_duallfp(hmm,emiss);

[state_durations] = compute_state_durations_seg(hmm_state_seq,Fs);

temp_p = zeros(2,length(dur_range));
temp_p(1,:) = geometric_pmf(1:length(dur_range),hmm.A(2,1));
temp_p(2,:) = geometric_pmf(1:length(dur_range),hmm.A(1,2));
hmm.p = temp_p;
for i = 1:hmm.K
    hmm.state(i).dur_type = 'hmm';
end

hmm.dur_range = dur_range;
hmm.max_state_dur = max_state_dur;
hmm.min_state_dur = min_state_dur;


%% for optimizing transition times
dsf_bb = 8;
Fs_bb = raw_Fs/dsf_bb;
bb_lcf = 0.05;
bb_hcf = 20;
hf_lcf = 20;
hf_hcf = 80;
hf_smooth = 0.05;
[obs_bb,t_axis,Fs_bb] = get_lf_features(raw_data,raw_Fs,Fs_bb,[bb_lcf bb_hcf]);
[obs_hf,t_axis,Fs_bb] = get_hf_features(raw_data,raw_Fs,Fs_bb,[hf_lcf hf_hcf],hf_smooth);
obs_bb = [obs_bb(:) obs_hf(:)];

pert_range_bb = [-0.15 0.15];

[bb_hmm_state_seq] = pert_optimize_hsmm_transitions_seg_dp_robust_duallfp... 
    (obs_bb,hmm,hmm_state_seq,gamma_hmm,Fs,Fs_bb,pert_range_bb);


