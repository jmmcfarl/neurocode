% this routine loads an example LFP recording (z-scored and resampled to 252Hz),
% and fits an HMM and a HSMM, displaying some diagnostic figures. 
% In this example, the UDS state sequences obtained from the HMM and HSMM are virtually
% identical.  

%specify path here
cur_dir = 'C:\WC_Germany\hsmm_uds_toolbox\hsmm_uds_toolbox';
cd(cur_dir)

%% Load test LFP data (Fs is sample-frequency =252Hz)
load ./test_lfp_data

%% Locate desynchronized epochs within the data (without clear UDS)
%% NOTE: you can assume that ALL data contains UDS by skipping this step and setting desynch_times = []; 
%set parameters
SO_thresh = -6; %threshold log relative power in the slow-oscillation band (0.2-1.5Hz)
HF_thresh = -2; %maximal power in the HF-band (>4Hz) to be considered a desynchronized epoch 
params.dsf = 1; %down-sample-factor for computing spectrogram
params.lcf = 0.1; %low-cutoff frequency 
params.hcf = 40; %high-cutoff frequency
params.movingwin = [15 2.5]; %windowing parameter for computing spectrogram (in seconds)

%compute desynchronized epochs and the spectrogram of the data
[desynch_times,desynch_ids,P,f,t] = hsmm_uds_locate_desynch_times(lfp,Fs,params,SO_thresh,HF_thresh);

%plot spectrogram and highlight desynchronized epochs
figure(1)
pcolor(t,f,log10(P)');shading flat, colorbar
ylim([.1 40])
set(gca,'yscale','log')
yl = ylim();
for i = 1:size(desynch_times,1)
   line([desynch_times(i,1) desynch_times(i,1)],yl,'color','w','linewidth',2);
   line([desynch_times(i,2) desynch_times(i,2)],yl,'color','w','linewidth',2);   
end
xlabel('Time','fontsize',14), ylabel('Frequency (Hz)','fontsize',14)

%% Compute signal features to use for UDS classification
hmm_dsf = 5; %down-sample-factor for HMM signal featurs
hmm_Fs = Fs/hmm_dsf; %sample-frequency of HMM signal features (Hz)
lf_cut_off_freqs = [0.05 2]; %for "low-frequency amplitude" features, these are the cutoff-frequencies for the bandpass filter
hf_cut_off_freqs = [20 80]; %for "high-frequency power" features, these are the cutoff-frequencies for the bandpass filter
smooth_sigma = 0.15; %smoothing sigma for computing HF power in seconds
% [sig_features,t_axis] = hsmm_uds_get_lf_features(lfp,Fs,hmm_Fs,lf_cut_off_freqs); %'low-frequency amplitude'
[sig_features,t_axis] = hsmm_uds_get_hf_features(lfp,Fs,hmm_Fs,hf_cut_off_freqs,smooth_sigma); %'high-frequency power'

%% LOCATE THE SEGMENTS CONTAINING UDS
T = length(sig_features); %number of samples
min_seg_dur = 60; %minimum duration of UDS segments used for analysis (in sec)
UDS_segs = hsmm_uds_get_uds_segments(desynch_times,hmm_Fs,T,min_seg_dur); %Nx2 matrix containing the index values of the beginning and end of each UDS segment

%% INITIALIZE THE HMM
params.meantype = 'variable'; %use 'variable' for time-varying state means. otherwise use 'fixed'
params.UDS_segs = UDS_segs;
params.movingwin = [5 1]; %moving window parameters [windowLength windowSlide](in seconds) for computing time-varying state means
[hmm] = hsmm_uds_initialize(sig_features,hmm_Fs,params);
% [hmm] = hsmm_uds_initialize(sig_features,hmm_Fs,params);

%% FIT AN HMM
hmm.min_mllik = -8; %set minimum value of the maximal log likelihood for either state before switching to robust estimator.  If you don't want to use this feature, set this to -Inf
hmm = hsmm_uds_train_hmm(hmm,sig_features);

%% DETERMINE THE VITERBI SEQUENCE
[state_seq,llik_best] = hsmm_uds_viterbi_hmm(hmm,sig_features);

%% OPTIMIZE STATE TRANSITIONS TO A BROADBAND SIGNAL WITH HIGHER SAMPLE-RATE
%compute broad-band signal
dsf_bb = 1; %down-sample factor (relative to original data)
bb_cf = [0.05 20]; %low and high cutoff frequencies for bandpass filtering (Hz)
pert_range_bb = [-0.15 0.15]; %range of allowed perturbations of state transitions (in sec)
[bb_features,t_axis,Fsd] = hsmm_uds_get_lf_features(lfp,Fs,Fs/dsf_bb,bb_cf); %'low-frequency amplitude'

%initialize state duration parameters
min_state_dur = 1; %minimum state duration (in samples)
max_state_dur = round(hmm_Fs*30); %maximum state duration in samples
dur_range = (1:max_state_dur)/hmm_Fs; %range of allowed state durations (in samples)
%set the state duration model for each hidden state to "HMM" ("Geometric")
for i = 1:hmm.K
    hmm.state(i).dur_type = 'hmm';
end
hmm.dur_range = dur_range;
hmm.max_state_dur = max_state_dur;
hmm.min_state_dur = min_state_dur;
[hmm_bb_state_seq] = hsmm_uds_pert_optimize_transitions(bb_features,hmm,state_seq,hmm.gamma,hmm.Fs,Fs/dsf_bb,pert_range_bb);

%create 0-1 state vector corresponding to the optimized state sequence
[new_seg_inds] = resample_uds_seg_inds(UDS_segs,hmm_Fs,Fs,length(lfp));
hmm_state_vec = nan(size(lfp));
for i = 1:size(new_seg_inds,1)
    hmm_state_vec(new_seg_inds(i,1):new_seg_inds(i,2)) = hmm_bb_state_seq{i};
end

%plot HMM state vector
figure(2)
plot(t_axis,lfp,'k'), hold on
plot(t_axis,hmm_state_vec,'r')
axis tight
xlabel('Time (s)','fontsize',14), ylabel('Amplitude (z)','fontsize',14)

%%
state_durations = hsmm_uds_compute_state_durations(state_seq,hmm_Fs);

for i = 1:hmm.K
    %estimate the empirical state duration pmf
    emp_pmf(i,:) = hist(state_durations{i},dur_range);
    emp_pmf(i,:) = emp_pmf(i,:)/sum(emp_pmf(i,:));  
    %estimate parameter of geometric pmf fit
    p_fit = 1/nanmedian(state_durations{i}*hmm_Fs);
    geo_pmf(i,:) = geometric_pmf(1:length(dur_range),p_fit);   
    [mu,lambda] = inverse_gauss_mlfit(state_durations{i});%compute the ML parameters of the inverse gaussian dist
    ig_pmf(i,:) = inverse_gaussian_pmf(dur_range,mu,lambda);%estimate the IG PMF    
    [alpha,beta] = gamma_mlfit(state_durations{i}); %compute the ML parameters of the gamma
    gam_pmf(i,:) = gamma_pmf(dur_range,alpha,beta); %estimate the gamma PMF
end
%plot empirical state duration distribution from HMM fit, along with
%paremetric model ML fits
figure(3)
subplot(2,1,1)
bar(dur_range,emp_pmf(1,:))
hold on
plot(dur_range,ig_pmf(1,:),'k',dur_range,gam_pmf(1,:),'r',dur_range,geo_pmf(i,:),'g','linewidth',2)
xlim([0 3.5])
title('DOWN State','fontsize',16)
xlabel('Duration (s)','fontsize',16)
ylabel('Relative frequency','fontsize',16)
legend('Empirical Distribution','Inverse Gaussian','Gamma','Geometric')
subplot(2,1,2)
bar(dur_range,emp_pmf(2,:))
hold on
plot(dur_range,ig_pmf(2,:),'k',dur_range,gam_pmf(2,:),'r',dur_range,geo_pmf(i,:),'g','linewidth',2)
xlim([0 3.5])
title('UP State','fontsize',16)
xlabel('Duration (s)','fontsize',16)
ylabel('Relative frequency','fontsize',16)

%%
%select duration model for each hidden state.  Choices are: 'inv_gauss', and 'gamma'
hmm.state(1).dur_type = 'inv_gauss';
hmm.state(2).dur_type = 'inv_gauss';
% hmm.state(1).dur_type = 'gamma';
% hmm.state(2).dur_type = 'gamma';
%initialize state duration model parameters
for i = 1:hmm.K
    if strcmp(hmm.state(i).dur_type,'inv_gauss')
        hmm.state(i).dur_pars = [mu lambda];
        hmm.P(i,:) = ig_pmf(i,:);
    else
        hmm.state(i).dur_pars = [alpha beta];
        hmm.P(i,:) = gam_pmf(i,:);
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

%% Fit HSMM
hsmm = hsmm_uds_train_hsmm(hmm,sig_features);

%% Compute HSMM Viterbi sequence
[hsmm_state_seq,max_lik] = hsmm_uds_viterbi_hsmm(hsmm,sig_features);

%% Optimize state transition times of the Viterbit sequence
[hsmm_bb_state_seq] = hsmm_uds_pert_optimize_transitions(bb_features,hsmm,hsmm_state_seq,hsmm.gamma,hsmm.Fs,Fs/dsf_bb,pert_range_bb);

%create 0-1 state vector corresponding to the optimized state sequence
hsmm_state_vec = nan(size(lfp));
for i = 1:size(new_seg_inds,1)
    hsmm_state_vec(new_seg_inds(i,1):new_seg_inds(i,2)) = hsmm_bb_state_seq{i};
end

%% add HSMM state sequence plot for comparison
figure(2)
plot(t_axis,hsmm_state_vec,'b')