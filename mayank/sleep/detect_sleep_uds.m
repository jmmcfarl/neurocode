clear all
cd ~/Analysis/Mayank/sleep/
load sleep_dirs
% load sleep_dirs_old

addpath(genpath('~/James_scripts/chronux/spectral_analysis/'))

%%
dd = 9;
cd(data(dd).dir)
load procData

if length(data(dd).ipsiLFPs) == 8
    poss_ctx_chs = 1:4
elseif length(data(dd).ipsiLFPs) == 16
    poss_ctx_chs = 1:8;
else
    error('No ipsi LFP chs');
end

% ipsi_csc = contra_csc;
%%

csc_Fs = 1/nanmedian(diff(csc_time));
mp_Fs = 1/nanmedian(diff(mp_t));

% mp_dsf = 5;
mp_dsf = 1;
csc_dsf = 1;
heka_dsf = 5;

%%
if ~isempty(heka_data)
heka_Fs = 1/nanmedian(diff(heka_time));
heka_data = decimate(heka_data,heka_dsf);
heka_Fsd = heka_Fs/heka_dsf;
[b,a] = butter(2,0.05/(heka_Fsd/2),'high');
heka_data = zscore(filtfilt(b,a,heka_data));
heka_time = downsample(heka_time,heka_dsf);
end

csc_Fsd = csc_Fs/csc_dsf;
[bb,aa] = butter(4,0.2/(csc_Fsd/2),'high');
for ii = 1:length(ipsi_csc)
    ipsi_csc{ii} = decimate(ipsi_csc{ii},csc_dsf);
    ipsi_csc_hp{ii} = filtfilt(bb,aa,ipsi_csc{ii});
    ipsi_csc{ii} = ipsi_csc{ii}/robust_std_dev(ipsi_csc_hp{ii});
    ipsi_csc_hp{ii} = ipsi_csc_hp{ii}/robust_std_dev(ipsi_csc_hp{ii});
end
csc_time = downsample(csc_time,csc_dsf);

mua_sm_sig = round(0.025*csc_Fsd);
ipsi_binned_spks = nan(length(csc_time),length(ipsi_mua_times));
ipsi_sm_rate = nan(length(csc_time),length(ipsi_mua_times));
for ii = 1:length(ipsi_mua_times)
    ipsi_binned_spks(:,ii) = hist(ipsi_mua_times{ii},csc_time);
    ipsi_sm_rate(:,ii) = zscore(jmm_smooth_1d_cor(ipsi_binned_spks(:,ii),mua_sm_sig));
end

%%
uu = find(~isnan(mp_t));
mp_interp = interp1(mp_t(uu),mp_d(uu),csc_time);
% uu = find(~isnan(heka_time));
% mp_interp = interp1(heka_time(uu),heka_data(uu),csc_time);
params.Fs = csc_Fsd;
params.tapers = [2 3];
params.fpass = [0 40];
movingwin = [20 10];
clear C

[b,a] = butter(2,[0.05]/(csc_Fs/2),'high');
mp_interp(isnan(mp_interp)) = 0;
% mp_interp = filtfilt(b,a,mp_interp);

use_cscs = 1:1:length(data(dd).ipsiLFPs);

clear C S2 S1
for cc = 1:length(use_cscs)
    cur_lfp = filtfilt(b,a,ipsi_csc{use_cscs(cc)});
%     cur_lfp = ipsi_csc{use_cscs(cc)};
[C{cc},phi,S12,S1,S2{cc},t,f]=cohgramc(mp_interp(:),cur_lfp(:),movingwin,params);
end

uds_range = [0.4 3];
hf_range = [20 80];
lf_range = [0 0.25];
use_range = [2 50];

uds_freqs = find(f >= uds_range(1) & f <= uds_range(2));
hf_freqs = find(f >= hf_range(1) & f <= hf_range(2));
lf_freqs = find(f >= lf_range(1) & f <= lf_range(2));
use_freqs = find(f >= use_range(1) & f <= use_range(2));

% mp_uds_pow = trapz(f(uds_freqs),log10(S1(:,uds_freqs)),2);
% mp_hf_pow = trapz(f(hf_freqs),log10(S1(:,hf_freqs)),2);
% mp_lf_pow = trapz(f(lf_freqs),log10(S1(:,lf_freqs)),2);
% mp_tot_pow = trapz(f(use_freqs),log10(S1(:,use_freqs)),2);
mp_uds_pow = log10(trapz(f(uds_freqs),(S1(:,uds_freqs)),2));
mp_hf_pow = log10(trapz(f(hf_freqs),(S1(:,hf_freqs)),2));
mp_lf_pow = log10(trapz(f(lf_freqs),(S1(:,lf_freqs)),2));
mp_tot_pow = log10(trapz(f(use_freqs),(S1(:,use_freqs)),2));

mp_tot_pow = (mp_tot_pow - median(mp_tot_pow))/robust_std_dev(mp_tot_pow);

clear lfp*pow
for cc = 1:length(use_cscs)
%    lfp_uds_pow(cc,:) = trapz(f(uds_freqs),log10(S2{cc}(:,uds_freqs)),2); 
%    lfp_hf_pow(cc,:) = trapz(f(hf_freqs),log10(S2{cc}(:,hf_freqs)),2); 
%    lfp_lf_pow(cc,:) = trapz(f(lf_freqs),log10(S2{cc}(:,lf_freqs)),2); 
   lfp_uds_pow(cc,:) = log10(trapz(f(uds_freqs),(S2{cc}(:,uds_freqs)),2)); 
   lfp_hf_pow(cc,:) = log10(trapz(f(hf_freqs),(S2{cc}(:,hf_freqs)),2)); 
   lfp_lf_pow(cc,:) = log10(trapz(f(lf_freqs),(S2{cc}(:,lf_freqs)),2)); 
end

%%
avg_uds_pow = mean(lfp_uds_pow,2);
avg_lf_pow = mean(lfp_lf_pow,2);
[~,ctx_ch] = max(avg_uds_pow(poss_ctx_chs)-0*avg_lf_pow(poss_ctx_chs));
ctx_ch = poss_ctx_chs(ctx_ch);

% lfp_uds_thresh = -2.5;
lfp_uds_thresh = -0.5;
% lfp_lf_max = 0.5;
lfp_lf_max = 1;
uds_epochs = (lfp_uds_pow(ctx_ch,:) >= lfp_uds_thresh);
lf_epochs = (lfp_lf_pow(ctx_ch,:) >= lfp_lf_max);

% has_uds = true(length(csc_time),1);
% for ii = 1:length(t)
%    if ~ismember(ii,uds_epochs)
%       curset = find(csc_time >= t(ii)-movingwin(1)/2 & csc_time <= t(ii) + movingwin(1)/2);
%       has_uds(curset) = false;
%    end
% end

%%
% f1 = figure();
% subplot(3,1,1); hold on
% imagesc(t,1:length(use_cscs),lfp_uds_pow); colorbar;
% plot(t,uds_epochs*4,'w','linewidth',2);
% line([0 t(end)],ctx_ch + [0 0],'color','k');
% axis tight
% xlim([0 t(end)]);
% 
% subplot(3,1,2); hold on
% imagesc(t,1:length(use_cscs),lfp_lf_pow); colorbar; 
% plot(t,lf_epochs*4,'w','linewidth',2);
% line([0 t(end)],ctx_ch + [0 0],'color','k');
% axis tight
% xlim([0 t(end)]);
% 
% subplot(3,1,3); hold on
% imagesc(t,f,log10(S2{ctx_ch})'); colorbar;
% plot(t,uds_epochs*4,'w','linewidth',2);
% ylim([0 5]);
% xlim([0 t(end)]);

%%
addpath('~/James_scripts/hsmm_uds_toolbox/')
removal_window = 0; %buffer window around desynch epochs

%mark points where the SO power crosses below threshold
desynch_indicator = ~uds_epochs | lf_epochs;
desynch_start_ids = 1+find(desynch_indicator(1:end-1) == 0 & desynch_indicator(2:end) == 1);
desynch_stop_ids = 1+find(desynch_indicator(1:end-1) == 1 & desynch_indicator(2:end) == 0);

%make sure you start and stop in desynchronized epochs in correct order
if ~isempty(desynch_start_ids)
    if isempty(desynch_stop_ids)
        desynch_stop_ids = length(t);
    else
        if desynch_start_ids(1) > desynch_stop_ids(1)
            desynch_start_ids = [1 desynch_start_ids];
        end
        if desynch_start_ids(end) > desynch_stop_ids(end)
            desynch_stop_ids = [desynch_stop_ids length(t)];
        end
    end
end

if length(desynch_start_ids) ~= length(desynch_stop_ids)
    disp('error start and stop IDs not equal!')
end

%compute the duration of each putative desynchronized epoch (in seconds)
desynch_durs = (desynch_stop_ids-desynch_start_ids)*movingwin(2);

%now make a window around desynchronized times for data exclusion
for w = 1:length(desynch_start_ids)
    if desynch_start_ids(w) <= removal_window
        desynch_start_ids(w) = 1;
    else
        desynch_start_ids(w) = desynch_start_ids(w)-removal_window;
    end
    if length(t)-desynch_stop_ids(w) <= removal_window
        desynch_stop_ids(w) = length(t);
    else
        desynch_stop_ids(w) = desynch_stop_ids(w)+removal_window;
    end
end

%now make sure there are no overlapping windows
bad_desynch_start = [];
for w = 2:length(desynch_start_ids)
    if desynch_start_ids(w) < desynch_stop_ids(w-1)
        bad_desynch_start = [bad_desynch_start w];
    end
end
desynch_start_ids(bad_desynch_start) = [];
desynch_stop_ids(bad_desynch_start-1) = [];

desynch_start_times = t(desynch_start_ids);
desynch_stop_times = t(desynch_stop_ids);
desynch_times = [desynch_start_times(:) desynch_stop_times(:)];
desynch_ids = round(desynch_times*csc_Fsd);

%% Compute signal features to use for UDS classification
hmm_dsf = 5; %down-sample-factor for HMM signal featurs
hmm_Fs = csc_Fs/hmm_dsf; %sample-frequency of HMM signal features (Hz)
lf_cut_off_freqs = [0.05 4]; %for "low-frequency amplitude" features, these are the cutoff-frequencies for the bandpass filter
[mp_features,t_axis] = hsmm_uds_get_lf_features(mp_interp,csc_Fs,hmm_Fs,lf_cut_off_freqs); %'low-frequency amplitude'
mp_features = mp_features/robust_std_dev(mp_features);

lf_cut_off_freqs = [0.1 2]; %for "low-frequency amplitude" features, these are the cutoff-frequencies for the bandpass filter
[lfp_features,t_axis] = hsmm_uds_get_lf_features(ipsi_csc{ctx_ch},csc_Fs,hmm_Fs,lf_cut_off_freqs); %'low-frequency amplitude'
% hf_cut_off_freqs = [20 80]; %for "high-frequency power" features, these are the cutoff-frequencies for the bandpass filter
% smooth_sigma = 0.15; %smoothing sigma for computing HF power in seconds
% [lfp_features,t_axis] = hsmm_uds_get_hf_features(ipsi_csc{ctx_ch},csc_Fs,hmm_Fs,hf_cut_off_freqs,smooth_sigma); %'high-frequency power'
lfp_features = lfp_features/robust_std_dev(lfp_features);
% lfp_features(lfp_features > 10) = 10;
%% LOCATE THE SEGMENTS CONTAINING UDS
T = length(mp_features); %number of samples
min_seg_dur = 20; %minimum duration of UDS segments used for analysis (in sec)
UDS_segs = hsmm_uds_get_uds_segments(desynch_times,hmm_Fs,T,min_seg_dur); %Nx2 matrix containing the index values of the beginning and end of each UDS segment

%% INITIALIZE THE HMM
params.meantype = 'variable'; %use 'variable' for time-varying state means. otherwise use 'fixed'
params.UDS_segs = UDS_segs;
params.movingwin = [10 1]; %moving window parameters [windowLength windowSlide](in seconds) for computing time-varying state means
[hmm_mp] = hsmm_uds_initialize(mp_features,hmm_Fs,params);
% [hmm_lfp] = hsmm_uds_initialize(lfp_features,hmm_Fs,params);

%% FIT AN HMM
hmm_mp.min_mllik = -8; %set minimum value of the maximal log likelihood for either state before switching to robust estimator.  If you don't want to use this feature, set this to -Inf
% hmm_lfp.min_mllik = -8; %set minimum value of the maximal log likelihood for either state before switching to robust estimator.  If you don't want to use this feature, set this to -Inf
fprintf('Fitting MP HMM\n');
hmm_mp = hsmm_uds_train_hmm(hmm_mp,mp_features);
fprintf('Fitting LFP HMM\n');
% hmm_lfp = hsmm_uds_train_hmm(hmm_lfp,lfp_features);

%% DETERMINE THE VITERBI SEQUENCE
[mp_state_seq,llik_best] = hsmm_uds_viterbi_hmm(hmm_mp,mp_features);
% [lfp_state_seq,llik_best] = hsmm_uds_viterbi_hmm(hmm_lfp,lfp_features);

%% OPTIMIZE STATE TRANSITIONS TO A BROADBAND SIGNAL WITH HIGHER SAMPLE-RATE
%compute broad-band signal
dsf_bb = 1; %down-sample factor (relative to original data)
pert_range_bb = [-0.15 0.15]; %range of allowed perturbations of state transitions (in sec)
bb_cf = [0.01 20]; %low and high cutoff frequencies for bandpass filtering (Hz)
[mp_bb_features,t_axis,Fsd] = hsmm_uds_get_lf_features(mp_interp,csc_Fs,csc_Fs/dsf_bb,bb_cf); %'low-frequency amplitude'
mp_bb_features = mp_bb_features/robust_std_dev(mp_bb_features);

bb_cf = [0.1 20]; %low and high cutoff frequencies for bandpass filtering (Hz)
[lfp_bb_features,t_axis,Fsd] = hsmm_uds_get_lf_features(ipsi_csc{ctx_ch},csc_Fs,csc_Fs/dsf_bb,bb_cf); %'low-frequency amplitude'
% bb_cf = [20 80]; %low and high cutoff frequencies for bandpass filtering (Hz)
% bb_smooth = 0.05;
% [lfp_bb_features,t_axis,Fsd] = hsmm_uds_get_hf_features(ipsi_csc{ctx_ch},csc_Fs,csc_Fs/dsf_bb,bb_cf,bb_smooth); %'low-frequency amplitude'
lfp_bb_features = lfp_bb_features/robust_std_dev(lfp_bb_features);
% lfp_bb_features(lfp_bb_features > 10) = 10;

%initialize state duration parameters
min_state_dur = 1; %minimum state duration (in samples)
max_state_dur = round(hmm_Fs*30); %maximum state duration in samples
dur_range = (1:max_state_dur)/hmm_Fs; %range of allowed state durations (in samples)
%set the state duration model for each hidden state to "HMM" ("Geometric")
for i = 1:hmm_mp.K
    hmm_mp.state(i).dur_type = 'hmm';
%     hmm_lfp.state(i).dur_type = 'hmm';
end
hmm_mp.dur_range = dur_range;
hmm_mp.max_state_dur = max_state_dur;
hmm_mp.min_state_dur = min_state_dur;
% hmm_lfp.dur_range = dur_range;
% hmm_lfp.max_state_dur = max_state_dur;
% hmm_lfp.min_state_dur = min_state_dur;
[mp_bb_state_seq] = hsmm_uds_pert_optimize_transitions(mp_bb_features,hmm_mp,mp_state_seq,hmm_mp.gamma,hmm_mp.Fs,csc_Fs/dsf_bb,pert_range_bb);
% [lfp_bb_state_seq] = hsmm_uds_pert_optimize_transitions(lfp_bb_features,hmm_lfp,lfp_state_seq,hmm_lfp.gamma,hmm_lfp.Fs,csc_Fs/dsf_bb,pert_range_bb);

%create 0-1 state vector corresponding to the optimized state sequence
[new_seg_inds] = resample_uds_seg_inds(UDS_segs,hmm_Fs,csc_Fs,length(mp_interp));
mp_state_vec = nan(size(mp_interp));
% lfp_state_vec = nan(size(mp_interp));
for i = 1:size(new_seg_inds,1)
    mp_state_vec(new_seg_inds(i,1):new_seg_inds(i,2)) = mp_bb_state_seq{i};
%     lfp_state_vec(new_seg_inds(i,1):new_seg_inds(i,2)) = lfp_bb_state_seq{i};
end

%%
% xl = [1e3 2.8e3];
% yl = [0 4];
% for cc = 1:length(use_cscs)
%     subplot(3,1,1)
%     pcolor(t,f,log(abs(S1))');shading flat
%     caxis([-5 0])
%     ylim(yl);
% %     xlim(xl);
% % set(gca,'yscale','log')
% 
%     subplot(3,1,2)
% pcolor(t,f,log(abs(S2{cc}))');shading flat
% % caxis([-6 1])
% caxis([25 32])
% % set(gca,'yscale','log')
%     ylim(yl);
% %     xlim(xl);
% 
%     subplot(3,1,3)
% pcolor(t,f,C{cc}');shading flat
% % set(gca,'yscale','log')
%     ylim(yl);
% %     xlim(xl);
% 
%     pause
%     clf
% end
