function [hmm] = parietal_get_hmm_lf_test_notrans(raw_data,raw_Fs,hcf,desynch_times,old_state_seq)
 
drive_letter = 'G';
addpath(strcat(drive_letter,':\Code\smoothing\software'))
addpath(strcat(drive_letter,':\Code\FullBNT-1.0.4\KPMstats\'))
addpath(strcat(drive_letter,':\Code\FullBNT-1.0.4\netlab3.3'))
addpath(strcat(drive_letter,':\WC_Germany\new_stellate_analysis\'))
addpath(strcat(drive_letter,':\WC_Germany\hsmm_state_detection'))
 
Fs_desired = 50; %(in Hz)
dsf = round(raw_Fs/Fs_desired); 
Fs = raw_Fs/dsf;
niqf = raw_Fs/2;
lf_cutoff_freqs = [0.05 hcf]; 
Fs_bb = raw_Fs/8;

windowSize = 50; %window duration to estimate time varying state means (s)
windowSlide = 1; %window spacing to estimate time varying state means (s)
num_states = 2; %number of hidden states (must be 2 in current implementation)
% state_mean_hcf = 0.01; %hcf for estimating time varying state means (0 for constant mean)
min_seg_dur = 60;

min_state_dur = 1;
max_state_dur = round(Fs*500);
min_seg_dur = 60;
dur_range = (1:max_state_dur)/Fs;


%% Get signal features
[obs_lf,t_axis,Fs] = get_lf_features(raw_data,raw_Fs,Fs_desired,lf_cutoff_freqs);
% [obs_lfd,t_axis,Fs] = get_lf_features(raw_data,raw_Fs,Fs_desired,[0.05 1]);
[obs_bb,t_axis,Fs_bb] = get_lf_features(raw_data,raw_Fs,Fs_bb,[0.05 40]);
T = length(obs_lf); 
delta = round(Fs*0.15/2); %number of bins on either side of each state transition to get rid of

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

%%
UDS_segs_c = UDS_segs;
for n = 1:length(old_state_seq)
    cur_state_seq = downsample(old_state_seq{n},5);
    cur_utrans = find(cur_state_seq(1:end-1) == 1 & cur_state_seq(2:end) == 2);
    cur_dtrans = find(cur_state_seq(1:end-1) == 2 & cur_state_seq(2:end) == 1);
    for j = 1:length(cur_utrans)
       if cur_utrans(j) > delta & cur_utrans(j) < length(cur_state_seq)-delta
           obs_lf(UDS_segs(n,1)+cur_utrans(j)-delta:UDS_segs(n,1)+cur_utrans(j)+delta) = nan;
       elseif cur_utrans(j) > delta
           obs_lf(UDS_segs(n,1)+cur_utrans(j)-delta:UDS_segs(n,2)) = nan;
       else
           obs_lf(UDS_segs(n,1):UDS_segs(n,1)+cur_utrans(j)+delta) = nan;  
       end
    end
    for j = 1:length(cur_dtrans)
       if cur_dtrans(j) > delta & cur_dtrans(j) < length(cur_state_seq)-delta
           obs_lf(UDS_segs(n,1)+cur_dtrans(j)-delta:UDS_segs(n,1)+cur_dtrans(j)+delta) = nan;
       elseif cur_dtrans(j) > delta
           obs_lf(UDS_segs(n,1)+cur_dtrans(j)-delta:UDS_segs(n,2)) = nan;
       else
           obs_lf(UDS_segs(n,1):UDS_segs(n,1)+cur_dtrans(j)+delta) = nan;  
       end
    end
    tot_chop = (length(cur_utrans)+length(cur_dtrans))*(2*delta+1);
    UDS_segs_c(n,2) = UDS_segs(n,2)-tot_chop;
    if n < size(UDS_segs_c,1)
        UDS_segs_c(n+1,:) = UDS_segs_c(n+1,:) - tot_chop;
    end
end
UDS_segs = UDS_segs_c;
%% initialize hsmm
meanparams.windowSize = windowSize;
meanparams.windowSlide = windowSlide;
% meanparams.state_mean_hcf = state_mean_hcf;
meanparams.meantype = 'variable';

% emiss = obs_lfd;
emiss = obs_lf;
emiss(isnan(emiss)) = [];
[hmm] = initialize_uds_hsmm_seg(emiss,Fs,meanparams,UDS_segs,drive_letter);
 
hmm.min_mllik = -2.;

[hmm,gamma] = hmm_uds_train_seg(hmm,emiss);
 

