function [hsmm_state_seq,hmm,state_durations] = parietal_get_hmm_lf_test(raw_data,raw_Fs,hcf,desynch_times)
 
drive_letter = 'G';
% addpath(strcat(drive_letter,':\Code\smoothing\software'))
% addpath(strcat(drive_letter,':\Code\FullBNT-1.0.4\KPMstats\'))
% addpath(strcat(drive_letter,':\Code\FullBNT-1.0.4\netlab3.3'))
% addpath(strcat(drive_letter,':\WC_Germany\new_stellate_analysis\'))
% addpath(strcat(drive_letter,':\WC_Germany\hsmm_state_detection'))

Fs_desired = 50; %(in Hz) this is the sampling frequency we want to work with
dsf = round(raw_Fs/Fs_desired); 
Fs = raw_Fs/dsf;
niqf = raw_Fs/2;
lf_cutoff_freqs = [0.05 hcf]; %use these as the bandpass filter cutoff frequencies
Fs_bb = raw_Fs/8; %for 'broadband' resampling of state sequence, use this sample frequency

windowSize = 50; %window duration to estimate time varying state means (s)
windowSlide = 10; %window spacing to estimate time varying state means (s)
num_states = 2; %number of hidden states (must be 2 in current implementation)

min_state_dur = 1; %minimum allowed state duration (in samples)
max_state_dur = round(Fs*500); %maximum allowed state duration in samples)
min_seg_dur = 60; %minimum UDS segment duration (in seconds)
dur_range = (1:max_state_dur)/Fs; %range of allowed state durations (in seconds)


%% Get signal features
[obs_lf,~,Fs] = get_lf_features(raw_data,raw_Fs,Fs_desired,lf_cutoff_freqs); %these are the low-frequency features used to classify HMM
[obs_bb,~,Fs_bb] = get_lf_features(raw_data,raw_Fs,Fs_bb,[0.05 40]); %these are the broad-band features used to realign state sequence
T = length(obs_lf); %number of samples in observation sequence 

%% compute UDS segments
UDS_segs = hsmm_uds_get_uds_segments(desynch_times,Fs,length(obs_lf),min_seg_dur);

%% initialize hsmm
% meanparams.windowSize = windowSize;
% meanparams.windowSlide = windowSlide;
% meanparams.meantype = 'variable';
params.meantype = 'variable';
params.UDS_segs = UDS_segs;
params.movingwin = [50 10];

emiss = obs_lf; %set emissions as the low-frequency amplitude
%initialize the HMM
% [hmm] = initialize_uds_hsmm_seg(emiss,Fs,meanparams,UDS_segs,drive_letter);
[hmm] = hsmm_uds_initialize(emiss,Fs,params);
 
% hmm.min_mllik = -2.; %this is the minimum max-likelihood before using the robust likelihood model
hmm.min_mllik = -8.; %this is the minimum max-likelihood before using the robust likelihood model
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

