function [hsmm_state_seq,hmm,state_durations] = parietal_get_hmm_hf_test(raw_data,raw_Fs,hf_smooth,desynch_times)
 
 %modified from parietal_get_hmm_lf_test
 
drive_letter = 'G';
% addpath(strcat(drive_letter,':\Code\smoothing\software'));
% addpath(strcat(drive_letter,':\Code\FullBNT-1.0.4\KPMstats\'));
% addpath(strcat(drive_letter,':\Code\FullBNT-1.0.4\netlab3.3'));
% addpath(strcat(drive_letter,':\WC_Germany\new_stellate_analysis\'));
% addpath(strcat(drive_letter,':\WC_Germany\hsmm_state_detection'));
 
Fs_desired = 50; %(in Hz)
dsf = round(raw_Fs/Fs_desired); 
Fs = raw_Fs/dsf;
niqf = raw_Fs/2;
hf_lcf = 20; %low cut-off for the low-freq filter
hf_hcf = 80; %high cut-off for the low-freq filter
hf_smooth2 = 0.05;
Fs_bb = raw_Fs/8;

windowSize = 50; %window duration to estimate time varying state means (s)
windowSlide = 1; %window spacing to estimate time varying state means (s)
num_states = 2; %number of hidden states (must be 2 in current implementation)
% state_mean_hcf = 0.01; %hcf for estimating time varying state means (0 for constant mean)
 min_seg_dur = 60;

min_state_dur = 1;
max_state_dur = round(Fs*500);
dur_range = (1:max_state_dur)/Fs;


%% Get signal features
[obs_hf,t_axis,Fs] = get_hf_features(raw_data,raw_Fs,Fs_desired,[hf_lcf hf_hcf],hf_smooth);
[obs_bb,t_axis,Fs_bb] = get_hf_features(raw_data,raw_Fs,Fs_bb,[hf_lcf hf_hcf],hf_smooth2);

T = length(obs_hf); 

%% compute UDS segments
UDS_segs = hsmm_uds_get_uds_segments(desynch_times,Fs,length(obs_hf),min_seg_dur);

%% initialize hsmm
% meanparams.windowSize = windowSize;
% meanparams.windowSlide = windowSlide;
% meanparams.meantype = 'variable';
params.meantype = 'variable';
params.UDS_segs = UDS_segs;
params.movingwin = [50 10];

emiss = obs_hf;
% [hmm] = initialize_uds_hsmm_seg(emiss,Fs,meanparams,UDS_segs,drive_letter);
[hmm] = hsmm_uds_initialize(emiss,Fs,params);

% hmm.min_mllik = -1.5;
hmm.min_mllik = -8;
hmm.windowSize = 50;

% [hmm,gamma] = hmm_uds_train_seg(hmm,emiss);
hmm = hsmm_uds_train_hmm(hmm,emiss);
 
% [hmm_state_seq] = hmm_uds_viterbi_seg(hmm,emiss);
[hmm_state_seq,llik_best] = hsmm_uds_viterbi_hmm(hmm,emiss);
 
% [state_durations] = compute_state_durations_seg(hmm_state_seq,Fs);
state_durations = hsmm_uds_compute_state_durations(hmm_state_seq,Fs);

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

% up_max_pert = 0.15;
% down_max_pert = 0.15;
pert_range_bb = [-0.3 0.3];
% [hsmm_state_seq] = pert_optimize_hsmm_transitions_seg_dp_robust...
%     (obs_bb,hmm,hmm_state_seq,gamma,Fs,Fs_bb,pert_range_bb);
% [hsmm_state_seq] = hsmm_pert_optimize_transitions(obs_bb,hmm,hmm_state_seq,hmm.gamma,hmm.Fs,Fs_bb,pert_range_bb);
[hsmm_state_seq] = hsmm_uds_pert_optimize_transitions_v2(obs_bb,hmm,hmm_state_seq,hmm.gamma,hmm.Fs,Fs_bb,pert_range_bb);

% [state_durations] = compute_state_durations_seg(hsmm_state_seq,Fs_bb);
state_durations = hsmm_uds_compute_state_durations(hsmm_state_seq,Fs_bb);

