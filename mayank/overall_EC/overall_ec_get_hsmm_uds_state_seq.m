function [hmm_bbstate_seq,hmm,Fs_bb,Fs,gamma] = ...
    overall_ec_get_hsmm_uds_state_seq(raw_data,raw_Fs,desynch_times,f_names)

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
num_states = 2; %number of hidden states (must be 2 in current implementation)

meanparams.windowSize = 70; %window duration to estimate time varying state means (s) (default 60
meanparams.windowSlide = 1; %window spacing to estimate time varying state means (s)
% meanparams.state_mean_hcf = 0.01; %hcf for estimating time varying state means (0 for constant mean)
meanparams.meantype = 'variable';

min_state_dur = 1;
max_state_dur = round(Fs*50); %default 60
min_seg_dur = 50; %default 60
dur_range = (1:max_state_dur)/Fs;
% dur_hist_range = linspace(0,10,100);

%% initializations
[obs_lf,t_axis,Fs] = get_lf_features(raw_data,raw_Fs,Fs_desired,[lf_lcf lf_hcf]);
[obs_lfd,t_axis,Fs] = get_lf_features(raw_data,raw_Fs,Fs_desired,[lf_lcf lf_hcf_d]);

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



%% initialize hmm
emiss = obs_lfd;
[hmm] = initialize_uds_hsmm_seg(emiss,Fs,meanparams,UDS_segs,drive_letter);
emiss = obs_lf;

hmm.min_mllik = -2;

[hmm, gamma_hmm] = hmm_uds_train_seg(hmm,emiss);

[hmm_state_seq,llik_best] = hmm_uds_viterbi_seg(hmm,emiss);

%% compute hmm state durations
smoothed_seq = thresh_state_smooth_seg(hmm_state_seq,Fs,100,100);
[state_durations] = compute_state_durations_seg(smoothed_seq,Fs);

for i = 1:hmm.K
    emp_pmf = hist(state_durations{i},dur_range);
    emp_pmf = emp_pmf/sum(emp_pmf);
    [mu,lambda] = inverse_gauss_mlfit(state_durations{i});
    ig_pmf = inverse_gaussian_pmf(dur_range,mu,lambda);
    [alpha,beta] = gamma_mlfit(state_durations{i});
    gam_pmf = gamma_pmf(dur_range,alpha,beta);
    ks_ig = max(abs(cumsum(emp_pmf) - cumsum(ig_pmf)));
    ks_gam = max(abs(cumsum(emp_pmf) - cumsum(gam_pmf)));
    if ks_ig > ks_gam
        hmm.state(i).dur_type = 'gamma';
        hmm.state(i).dur_pars = [alpha beta];
        hmm.P(i,:) = gam_pmf;
    else
        hmm.state(i).dur_type = 'inv_gauss';
        hmm.state(i).dur_pars = [mu lambda];
        hmm.P(i,:) = ig_pmf;       
    end
end

hmm.dur_range = dur_range;
hmm.max_state_dur = max_state_dur;
hmm.min_state_dur = min_state_dur;

% enforce any minimum state duration and renormalize
for i = 1:hmm.K
    if min_state_dur > 1
        hmm.P(i,1:min_state_dur-1) = zeros(1,min_state_dur-1);
        hmm.P(i,:) = hmm.P(i,:)/sum(hmm.P(i,:));
    end
end

%% train hsmm
% [hmm,gamma,hmm_window_post]=hsmm_uds_train_seg(hmm,emiss);
% 
% [hsmm_state_seq,max_lik] = hsmm_uds_viterbi_seg(hmm,emiss);
% 
% hmm.posteriors = hmm_window_post;
% hmm.max_lik = max_lik;
% hmm.gamma = gamma;

% hsmm_state_seq = hmm_state_seq;
gamma = gamma_hmm;
%% optimizing transition times
dsf_bb = 8;
Fs_bb = raw_Fs/dsf_bb;
bb_lcf = 0.05;
bb_hcf = 20;
[obs_bb,t_axis,Fs_bb] = get_lf_features(raw_data,raw_Fs,Fs_bb,[bb_lcf bb_hcf]);

pert_range_bb = [-0.15 0.15];
[hmm_bbstate_seq] = pert_optimize_hsmm_transitions_seg_dp_robust... 
    (obs_bb,hmm,hmm_state_seq,gamma,Fs,Fs_bb,pert_range_bb);

%% plot meanfuns
% [obs_dist,obsrange,win_t] = compute_slide_dist(obs_lf,meanparams.windowSize,meanparams.windowSlide,Fs);
[obs_dist,obsrange,win_t] = compute_slide_dist(obs_lf,70,0.5,Fs);
offset = 0;

figure
imagesc(win_t,obsrange,log(obs_dist'))
caxis([-4 0])
hold on
for i = 1:hmm.Nsegs
    t_axis = (hmm.UDS_segs(i,1):hmm.UDS_segs(i,2))/Fs;
    plot(t_axis,hmm.state(1).meanfun{i},'w')
    plot(t_axis,hmm.state(2).meanfun{i},'w')
end
ylim([-2 2.5])
for i = 1:size(desynch_times,1)
    line([desynch_times(i,1) desynch_times(i,1)]+offset,[-2 2.5])
    line([desynch_times(i,2) desynch_times(i,2)]-offset,[-2 2.5])
end

t_names = ['G:\WC_Germany\overall_EC\meanfuns\' f_names];
print('-dpng',t_names), close


