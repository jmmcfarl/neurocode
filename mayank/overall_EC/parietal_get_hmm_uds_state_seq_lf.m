function [bb_hmm_state_seq,hmm_state_seq,hmm,Fs] = parietal_get_hmm_uds_state_seq_lf...
    (raw_data,raw_Fs,desynch_times,f_names)

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
num_states = 2; %number of hidden states (must be 2 in current implementation)

meanparams.windowSize = 50; %window duration to estimate time varying state means (s)
meanparams.windowSlide = 1; %window spacing to estimate time varying state means (s)
meanparams.meantype = 'variable';

min_state_dur = 1;
max_state_dur = round(Fs*30);
min_seg_dur = 60;
dur_range = (1:max_state_dur)/Fs;
dur_hist_range = linspace(0,10,100);

%% initializations
[obs_lf,t_axis,Fs] = get_lf_features(raw_data,raw_Fs,Fs_desired,[lf_lcf lf_hcf]);

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
emiss = obs_lf;
[hmm] = initialize_uds_hsmm_seg(emiss,Fs,meanparams,UDS_segs,drive_letter);
hmm.min_state_dur = min_state_dur;
hmm.max_state_dur = max_state_dur;
hmm.dur_range = dur_range;
hmm.min_mllik = -1.5; %set minimum value of the maximal log likelihood for either state before switching to robust estimator

[hmm, gamma_hmm] = hmm_uds_train_seg(hmm,emiss);

[hmm_state_seq,llik_best] = hmm_uds_viterbi_seg(hmm,emiss);

%%
dsf_bb = 8;
Fs_bb = raw_Fs/dsf_bb;
bb_lcf = 0.05;
bb_hcf = 20;
[obs_bb,t_axis,Fs_bb] = get_lf_features(raw_data,raw_Fs,Fs_bb,[bb_lcf bb_hcf]);

pert_range_bb = [-0.15 0.15];
% [bb_hmm_state_seq] = pert_optimize_hmm_transitions_seg_dp_robust... 
%     (obs_bb,hmm,hmm_state_seq,gamma_hmm,Fs,Fs_bb,pert_range_bb);

hmm_temp = hmm;
% temp_p = zeros(size(hmm.P));
temp_p(1,:) = geometric_pmf(1:length(dur_range),hmm_temp.A(2,1));
temp_p(2,:) = geometric_pmf(1:length(dur_range),hmm_temp.A(1,2));
hmm_temp.P = temp_p;
for i = 1:hmm.K
    hmm_temp.state(i).dur_type = 'hmm';
end

[bb_hmm_state_seq] = pert_optimize_hsmm_transitions_seg_dp_robust... 
    (obs_bb,hmm_temp,hmm_state_seq,gamma_hmm,Fs,Fs_bb,pert_range_bb);


%%
[obs_dist,obsrange,win_t] = compute_slide_dist(emiss,meanparams.windowSize,meanparams.windowSlide,Fs);
figure
pcolor(win_t,obsrange,log(obs_dist'));shading flat
caxis([-4 0])
hold on
for i = 1:hmm.Nsegs
    t_axis = (hmm.UDS_segs(i,1):hmm.UDS_segs(i,2))/hmm.Fs;
    plot(t_axis,hmm.state(1).meanfun{i},'w')
    plot(t_axis,hmm.state(2).meanfun{i},'w')
end
t_names = ['G:\WC_Germany\overall_EC\state_meanfuns\lf_' f_names];
print('-dpng',t_names), close

