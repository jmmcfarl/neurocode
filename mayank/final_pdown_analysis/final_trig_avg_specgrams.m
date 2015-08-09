
clear all
% addpath('C:/WC_Germany/parietal_cortical_2010/');
% addpath('C:/Code/general/');
% addpath('C:/Code/general_functions/');
% addpath('C:/WC_Germany/hsmm_state_detection/');
% addpath('C:/WC_Germany/persistent_downs/');
% addpath(genpath('C:/Code/figuremaker/'));
addpath('~/James_scripts/mayank/parietal_cortical_2010/');
addpath('~/James_scripts/mayank/persistent_9_27_2010/');
addpath('~/James_scripts/mayank/hsmm_state_detection/');
addpath('~/James_scripts/mayank/persistent_downs/');
addpath('~/James_scripts/mayank/final_pdown_analysis/');

% load C:/WC_Germany/final_pdown_analysis/compiled_data.mat
load ~/Analysis/Mayank/final_pdown_analysis/compiled_data.mat

% fig_dir = 'C:\WC_Germany\final_pdown_analysis\figures\';

min_rec_dur = 500; %minimum recording duration (in sec)
data_ids = [data(:).id];
used_dirs = find([data(:).ep] > min_rec_dur & [data(:).dp] > min_rec_dur); %only analyze recs where we have at least this much data to work with
% used_dirs([data(used_dirs).id] == 76) = []; %this rec does not have usable LFP UDS
used_dirs(ismember(data_ids(used_dirs),no_cell) | ismember(data_ids(used_dirs),unclear_uds)) = [];
% used_dirs(ismember(data_ids(used_dirs),unclear_uds)) = [];

data = data(used_dirs);


% load C:/WC_Germany/final_pdown_analysis/fin_pdown_core_analysis.mat
% load C:/WC_Germany/final_pdown_analysis/mua_classification2.mat
load ~/Analysis/Mayank/final_pdown_analysis/fin_pdown_core_analysis_fin.mat
load ~/Analysis/Mayank/final_pdown_analysis/mua_classification_fin.mat

% %include channels where hpc MUA peak appears on ch5
% peak_hpcmua_loc([27 29 67 68 79]) = 5;

peak_hpcmua_loc = peak_hpcmua_loc(used_dirs);
peak_hpcmua_rate = peak_hpcmua_rate(used_dirs);

if length(core_data) ~= length(data)
    error('Data mismatch');
end
% load C:/WC_Germany/final_pdown_analysis/final_cortical_state_data.mat

%%
raw_Fs = 2016;
dsf = 6;
usfac = 8/dsf; %ratio of original dsf to current one
Fsd = raw_Fs/dsf;
niqf = raw_Fs/2;
lcf = 0.05;hcf = 10; %bandpass filter for LF-amp
lcf_hf = 30; hcf_hf = 100; %bandpass filter for HF-power
hcf_sm = 0.025; %power smooth sigma
rate_sm = round(Fsd*0.05); %MUA rate smooth sigma

%define parameters of notch filters at 50 and 100Hz to remove line noise
wo = 50/(Fsd/2);
bw = wo/50;
[b_notch,a_notch] = iirnotch(wo,bw);
wo = 100/(Fsd/2);
[b_notch2,a_notch2] = iirnotch(wo,bw);

%lag windows for computing trig avgs
backlag = 1*Fsd;
forwardlag = 2*Fsd;

lags = (-backlag:forwardlag);
min_samps = 5; %minimum number of samples to compute density at a particular time point

%parameters of wavelet transform
n_freqs = 35; %number of scales
min_freq = 4; %min pseudo-freq
max_freq = 110; %max pseudo-freq
desired_wfreqs = linspace(min_freq,max_freq,n_freqs);
scales = Fsd./(desired_wfreqs);
wfreqs = scal2frq(scales,'cmor1-1',1/Fsd);

%%
for d = 1:length(data)
%     cd(data(d).dir)
    cd(data(d).new_dir)
    pwd
    spec_data(d).data_id = data(d).id;
    
    load ./used_data lf7 wcv_minus_spike

    %get synced timestamps
    load ./sync_times.mat
    synct_d = downsample(synct,dsf);
    
    %     if data(d).hpc_lfp == 3
    load ./used_data lf2 lf3 lf4 lf5
    if data(d).is_old_type
        lf3 = lf3 + lf5; %redefine LF3 wrt gnd
    end
    hpc_lfp = lf3; %default is LF3 is cell-layer hpc LFP
    dhpc_lfp = lf2; %and LF2 is deep hpc LFP
    %     end
    
    %if peak of MUA profile was not at the default position, adjust which
    %channels were pulling hpc LFP signals from
    if peak_hpcmua_loc(d) == 5
        hpc_lfp = lf5;
        dhpc_lfp = lf4;
    elseif peak_hpcmua_loc(d) == 4
        hpc_lfp = lf4;
        dhpc_lfp = lf3;
    end

    [lfp_lf,t_axis] = get_lf_features(lf7,raw_Fs,Fsd,[lcf hcf]);

    %down-sample hpc LFP, then notch filter, then zscore
    hpc_ds = decimate(hpc_lfp,dsf);
    hpc_ds = filtfilt(b_notch,a_notch,hpc_ds);
    hpc_ds = filtfilt(b_notch2,a_notch2,hpc_ds);
    hpc_ds = zscore(hpc_ds);
    
    %ditto for ctx LFP
    lfp_ds = decimate(lf7,dsf);
    lfp_ds = filtfilt(b_notch,a_notch,lfp_ds);
    lfp_ds = filtfilt(b_notch2,a_notch2,lfp_ds);
    lfp_ds = zscore(lfp_ds);
    
    %ditto for MP
    wcv_ds = decimate(wcv_minus_spike,dsf);
    wcv_ds = filtfilt(b_notch,a_notch,wcv_ds);
    wcv_ds = filtfilt(b_notch2,a_notch2,wcv_ds);
    wcv_ds = zscore(wcv_ds);
    
    %chop data after last usabe sample
    end_time = min(data(d).ep,data(d).dp);
    ep = find(t_axis >= end_time,1);
    if ~isempty(ep)
        synct_d(ep+1:end) = []; lfp_lf(ep+1:end) = []; wcv_ds(ep+1:end) = []; t_axis(ep+1:end) = [];
        lfp_ds(ep+1:end) = []; hpc_ds(ep+1:end) = [];
    else
        ep = length(t_axis);
    end
    
        %% get hsmm state seq data
    load ./pa_hsmm_state_seq_combined_fin_nd.mat
    load ./pa_hsmm_state_seq7_combined_fin_nd.mat
    hsmm_ctx = hsmm7;
    lfp_state_seq = hsmm_bbstate_seq7;
    mp_state_seq = hsmm_bbstate_seq;
    %resample HMM state seqs onto our current time axis (usfac is the ratio
    %of original dsf to current dsf)
    for i = 1:length(lfp_state_seq)
        temp = ceil((1/usfac):(1/usfac):length(lfp_state_seq{i}));
        lfp_state_seq{i} = lfp_state_seq{i}(temp);
        mp_state_seq{i} = mp_state_seq{i}(temp);
    end

    [new_seg_inds] = resample_uds_seg_inds_v2(hsmm.UDS_segs,hsmm.Fs,Fsd,mp_state_seq);
    dur_uds = sum(diff(new_seg_inds,[],2))/Fsd; %total duration of UDS segs
    [up_trans_inds_mp,down_trans_inds_mp] = compute_state_transitions_seg(new_seg_inds,mp_state_seq);
    [up_trans_inds_lfp,down_trans_inds_lfp] = compute_state_transitions_seg(new_seg_inds,lfp_state_seq);
    
    bad_mp_states = find(up_trans_inds_mp > ep | down_trans_inds_mp > ep);
    up_trans_inds_mp(bad_mp_states) = []; down_trans_inds_mp(bad_mp_states) = [];
    bad_lfp_states = find(up_trans_inds_lfp > ep | down_trans_inds_lfp > ep);
    up_trans_inds_lfp(bad_lfp_states) = []; down_trans_inds_lfp(bad_lfp_states) = [];
    
    %% get LFP cortical UDS period vector
    load ./allEC_ctx_period_data_hsmm.mat
    lfp_period_vec = nan(size(wcv_ds));
    for i = 1:length(lf8_period_f) %resample onto new time axis
        temp = ceil((1/usfac):(1/usfac):length(lf8_period_f{i}));
        lf8_period_f{i} = lf8_period_f{i}(temp);
    end
    for i = 1:size(new_seg_inds,1)
        cur_inds = new_seg_inds(i,1):new_seg_inds(i,2);
        cur_inds_used = find(cur_inds <= ep);
        lfp_period_vec(cur_inds(cur_inds_used)) = lf8_period_f{i}(cur_inds_used);
    end
    
    %% get state vectors
    mp_state_number = nan(length(wcv_ds),1);
    mp_state_vec = nan(length(wcv_ds),1);
    lfp_state_number = nan(length(wcv_ds),1);
    lfp_state_vec = nan(length(wcv_ds),1);
    for i = 1:length(up_trans_inds_mp)-1
        mp_state_vec(up_trans_inds_mp(i):down_trans_inds_mp(i)) = 1;
        mp_state_vec(down_trans_inds_mp(i):up_trans_inds_mp(i+1)) = 0;
        mp_state_number(up_trans_inds_mp(i):down_trans_inds_mp(i)) = 2*(i-1)+1;
        mp_state_number(down_trans_inds_mp(i):up_trans_inds_mp(i+1)) = 2*(i-1) + 2;
    end
    
    for i = 1:length(up_trans_inds_lfp)-1
        lfp_state_vec(up_trans_inds_lfp(i):down_trans_inds_lfp(i)) = 1;
        lfp_state_vec(down_trans_inds_lfp(i):up_trans_inds_lfp(i+1)) = 0;
        lfp_state_number(up_trans_inds_lfp(i):down_trans_inds_lfp(i)) = 2*(i-1)+1;
        lfp_state_number(down_trans_inds_lfp(i):up_trans_inds_lfp(i+1)) = 2*(i-1) + 2;
    end
    mp_state_number(isnan(lfp_period_vec)) = nan;
    lfp_state_number(isnan(lfp_period_vec)) = nan;
    mp_state_vec(isnan(lfp_period_vec)) = nan;
    lfp_state_vec(isnan(lfp_period_vec)) = nan;
    
    %% compute corresponding state transitions and transition lags
    [corresp_lfp_upinds,corresp_lfp_downinds] = find_corresponding_state_transitions_lookback(...
        up_trans_inds_mp,down_trans_inds_mp,up_trans_inds_lfp,down_trans_inds_lfp);

    %% get wavelet spectrograms
    uds_inds = find(~isnan(mp_state_vec)); %set of inds during UDS segs
    
%     lfp_specgram = abs(cwt(lfp_ds,scales,'cmor1-1'))';
%     hpc_specgram = abs(cwt(hpc_ds,scales,'cmor1-1'))';
%     mp_specgram = abs(cwt(wcv_ds,scales,'cmor1-1'))';
    lfp_specgram = log(abs(cwt(lfp_ds,scales,'cmor1-1')).^2)'; %log power
    hpc_specgram = log(abs(cwt(hpc_ds,scales,'cmor1-1')).^2)';
    mp_specgram = log(abs(cwt(wcv_ds,scales,'cmor1-1')).^2)';

    %compute avg and SD of spectra across time
    spec_data(d).avg_lfp_spectrum = mean(lfp_specgram(uds_inds,:));
    spec_data(d).std_lfp_spectrum = std(lfp_specgram(uds_inds,:));
    spec_data(d).avg_hpc_spectrum = mean(hpc_specgram(uds_inds,:));
    spec_data(d).std_hpc_spectrum = std(hpc_specgram(uds_inds,:));
    spec_data(d).avg_mp_spectrum = mean(mp_specgram(uds_inds,:));
    spec_data(d).std_mp_spectrum = std(mp_specgram(uds_inds,:));
    
    %zscore normalize specgrams 
    lfp_specgram = bsxfun(@rdivide,bsxfun(@minus,lfp_specgram,spec_data(d).avg_lfp_spectrum),spec_data(d).std_lfp_spectrum);
    hpc_specgram = bsxfun(@rdivide,bsxfun(@minus,hpc_specgram,spec_data(d).avg_hpc_spectrum),spec_data(d).std_hpc_spectrum);
    mp_specgram = bsxfun(@rdivide,bsxfun(@minus,mp_specgram,spec_data(d).avg_mp_spectrum),spec_data(d).std_mp_spectrum);

%     z_thresh = 5;
%     lf8_specgram(lf8_specgram > z_thresh) = z_thresh;
%     lf8_specgram(lf8_specgram < -z_thresh) = -z_thresh;
%     lf8_specgram = zscore(lf8_specgram);
%     
%     hpc_specgram(hpc_specgram > z_thresh) = z_thresh;
%     hpc_specgram(hpc_specgram < -z_thresh) = -z_thresh;
%     hpc_specgram = zscore(hpc_specgram);

    
    %% FOR CORTICAL UP TRANSITIONS
    
    %get skipped and non-skipped ctx up trans, but get rid of up-trans that
    % 'could not have been' skipped by MP
    ctx_ups_mp_up = find(mp_state_vec(up_trans_inds_lfp) == 1); %ctx ups when mp is already up
    skipped_lfp_ups = [core_data(d).skipped_lfp_ups];
    non_skipped_lfp_ups = [core_data(d).non_skipped_lfp_ups];
    skipped_lfp_ups(ismember(skipped_lfp_ups,ctx_ups_mp_up)) = [];
    non_skipped_lfp_ups(ismember(non_skipped_lfp_ups,ctx_ups_mp_up)) = [];
    
    %trig avg specgrams for non-skipped ctx up-trans
    non_sk_utrig_mp_spec = get_event_trig_avg(mp_specgram,up_trans_inds_lfp(non_skipped_lfp_ups),backlag,forwardlag);
    non_sk_utrig_lfp_spec = get_event_trig_avg(lfp_specgram,up_trans_inds_lfp(non_skipped_lfp_ups),backlag,forwardlag);
    non_sk_utrig_hpc_spec = get_event_trig_avg(hpc_specgram,up_trans_inds_lfp(non_skipped_lfp_ups),backlag,forwardlag);
    if length(skipped_lfp_ups) > min_samps %if there are enough examples, do for skipped ctx ups
        sk_utrig_mp_spec = get_event_trig_avg(mp_specgram,up_trans_inds_lfp(skipped_lfp_ups),backlag,forwardlag);
        sk_utrig_lfp_spec = get_event_trig_avg(lfp_specgram,up_trans_inds_lfp(skipped_lfp_ups),backlag,forwardlag);
        sk_utrig_hpc_spec = get_event_trig_avg(hpc_specgram,up_trans_inds_lfp(skipped_lfp_ups),backlag,forwardlag);
    else
        sk_utrig_mp_spec = nan(length(lags),length(wfreqs));
        sk_utrig_lfp_spec = nan(length(lags),length(wfreqs));
        sk_utrig_hpc_spec = nan(length(lags),length(wfreqs));
    end
        
    spec_data(d).non_sk_utrig_mp_spec = non_sk_utrig_mp_spec;
    spec_data(d).non_sk_utrig_lfp_spec = non_sk_utrig_lfp_spec;
    spec_data(d).non_sk_utrig_hpc_spec = non_sk_utrig_hpc_spec;
    spec_data(d).sk_utrig_mp_spec = sk_utrig_mp_spec;
    spec_data(d).sk_utrig_lfp_spec = sk_utrig_lfp_spec;
    spec_data(d).sk_utrig_hpc_spec = sk_utrig_hpc_spec;
    
    %% FOR CORTICAL DOWN TRANSITIONS
    %get set of skipped and non-skipped ctx down-trans (excluding those
    %that could not have been skipped)
    ctx_downs_mp_down = find(mp_state_vec(down_trans_inds_lfp) == 0);
    skipped_lfp_downs = [core_data(d).skipped_lfp_downs];
    non_skipped_lfp_downs = [core_data(d).non_skipped_lfp_downs];
    skipped_lfp_downs(ismember(skipped_lfp_downs,ctx_downs_mp_down)) = [];
    non_skipped_lfp_downs(ismember(non_skipped_lfp_downs,ctx_downs_mp_down)) = [];
    
    %ctx non-skipped down-trig avg specgrams
    non_sk_dtrig_mp_spec = get_event_trig_avg(mp_specgram,down_trans_inds_lfp(non_skipped_lfp_downs),backlag,forwardlag);
    non_sk_dtrig_lfp_spec = get_event_trig_avg(lfp_specgram,down_trans_inds_lfp(non_skipped_lfp_downs),backlag,forwardlag);
    non_sk_dtrig_hpc_spec = get_event_trig_avg(hpc_specgram,down_trans_inds_lfp(non_skipped_lfp_downs),backlag,forwardlag);
    if length(skipped_lfp_downs) > min_samps %if there were enough examples, do for skipped downs
        sk_dtrig_mp_spec = get_event_trig_avg(mp_specgram,down_trans_inds_lfp(skipped_lfp_downs),backlag,forwardlag);
        sk_dtrig_lfp_spec = get_event_trig_avg(lfp_specgram,down_trans_inds_lfp(skipped_lfp_downs),backlag,forwardlag);
        sk_dtrig_hpc_spec = get_event_trig_avg(hpc_specgram,down_trans_inds_lfp(skipped_lfp_downs),backlag,forwardlag);
    else
        sk_dtrig_mp_spec = nan(length(lags),length(wfreqs));
        sk_dtrig_lfp_spec = nan(length(lags),length(wfreqs));
        sk_dtrig_hpc_spec = nan(length(lags),length(wfreqs));
    end
        
    spec_data(d).non_sk_dtrig_mp_spec = non_sk_dtrig_mp_spec;
    spec_data(d).non_sk_dtrig_lfp_spec = non_sk_dtrig_lfp_spec;
    spec_data(d).non_sk_dtrig_hpc_spec = non_sk_dtrig_hpc_spec;
    spec_data(d).sk_dtrig_mp_spec = sk_dtrig_mp_spec;
    spec_data(d).sk_dtrig_lfp_spec = sk_dtrig_lfp_spec;
    spec_data(d).sk_dtrig_hpc_spec = sk_dtrig_hpc_spec;
    
    %%
end

%%
% cd C:\WC_Germany\final_pdown_analysis\
% save final_trig_spec_fin_noconst_newpeaks_wch5 spec_data wfreqs lags Fsd dsf
% save final_trig_spec_log_differential spec_data wfreqs lags Fsd dsf
cd ~/Analysis/Mayank/final_pdown_analysis/
save final_trig_spec_fin_noconst_newpeaks_wch5_fin spec_data wfreqs lags Fsd dsf

%%
% % load final_trig_avg_data3 trig_data
% 
% data_ids = [data(:).id];
% l3mec = find(strcmp({data.loc},'MEC'));
% l3lec = find(strcmp({data.loc},'LEC'));
% l3mec(~ismember(data_ids(l3mec),clear_l3)) = [];
% 
% fract_rt2_ups = [core_data(:).fract_rt2_ups];
% fract_rt2_downs = [core_data(:).fract_rt2_downs];
% 
% % n_skipped_ups = [trig_data(:).n_skipped_ups];
% % n_skipped_downs = [trig_data(:).n_skipped_downs];
% n_pers_downs = cellfun(@(x) length(x),{core_data(:).rt2_downs});
% n_pers_ups = cellfun(@(x) length(x),{core_data(:).rt2_ups});
% 
% min_states = 5;
% % l3mec_pups = l3mec(n_skipped_downs(l3mec) >= min_states);
% % l3mec_pdowns = l3mec(n_skipped_ups(l3mec) >= min_states);
% l3mec_pups = l3mec(n_pers_ups(l3mec) >= min_states);
% l3mec_pdowns = l3mec(n_pers_downs(l3mec) >= min_states);
% 
% interp_f = linspace(min(wfreqs),max(wfreqs),100)';
% %%
% used_data = l3mec_pdowns;
% xr = [-0.5 1.5];
% 
% non_sk_utrig_lf8_spec = reshape([spec_data(used_data).non_sk_utrig_lfp_spec],[length(lags) length(wfreqs) length(used_data)]);
% sk_utrig_lf8_spec = reshape([spec_data(used_data).sk_utrig_lfp_spec],[length(lags) length(wfreqs) length(used_data)]);
% % dom_utrig_lf8_spec = (non_sk_utrig_lf8_spec - sk_utrig_lf8_spec);
% 
% % ufreqs = find(wfreqs >= 40);
% % non_sk_utrig_lf8_hf = squeeze(mean(non_sk_utrig_lf8_spec(:,ufreqs,:),2));
% % sk_utrig_lf8_hf = squeeze(mean(sk_utrig_lf8_spec(:,ufreqs,:),2));
% % non_sk_utrig_hpc_hf = squeeze(mean(non_sk_utrig_hpc_spec(:,ufreqs,:),2));
% % sk_utrig_hpc_hf = squeeze(mean(sk_utrig_hpc_spec(:,ufreqs,:),2));
% 
% nonsk_specgram = squeeze(nanmean(non_sk_utrig_lf8_spec,3))';
% interp_nonsk_specgram = interp1(wfreqs,nonsk_specgram,interp_f);
% sk_specgram = squeeze(nanmean(sk_utrig_lf8_spec,3))';
% interp_sk_specgram = interp1(wfreqs,sk_specgram,interp_f);
% 
% h1 = figure;
% subplot(2,1,1);
% % pcolor(lags/Fsd,wfreqs,squeeze(nanmean(non_sk_utrig_lf8_spec,3))');shading flat
% imagesc(lags/Fsd,wfreqs,interp_nonsk_specgram);set(gca,'ydir','normal');
% xlim(xr);
% xlabel('Time (s)');
% ylabel('Frequency (Hz)');
% title('Non-skipped');
% colorbar
% ca = caxis();
% subplot(2,1,2);
% % pcolor(lags/Fsd,wfreqs,squeeze(nanmean(sk_utrig_lf8_spec,3))');shading flat
% imagesc(lags/Fsd,wfreqs,interp_sk_specgram);set(gca,'ydir','normal');
% xlim(xr);
% caxis(ca);
% xlabel('Time (s)');
% ylabel('Frequency (Hz)');
% title('Skipped');
% colorbar
% 
% % hh1 = figure();
% % pcolor(lags/Fsd,wfreqs,squeeze(nanmean(dom_utrig_lf8_spec,3))');shading flat
% 
% 
% used_data = l3mec_pdowns;
% 
% non_sk_utrig_mp_spec = reshape([spec_data(used_data).non_sk_utrig_mp_spec],[length(lags) length(wfreqs) length(used_data)]);
% sk_utrig_mp_spec = reshape([spec_data(used_data).sk_utrig_mp_spec],[length(lags) length(wfreqs) length(used_data)]);
% 
% nonsk_specgram = squeeze(nanmean(non_sk_utrig_mp_spec,3))';
% interp_nonsk_specgram = interp1(wfreqs,nonsk_specgram,interp_f);
% sk_specgram = squeeze(nanmean(sk_utrig_mp_spec,3))';
% interp_sk_specgram = interp1(wfreqs,sk_specgram,interp_f);
% 
% h2=figure;
% subplot(2,1,1);
% % pcolor(lags/Fsd,wfreqs,squeeze(nanmean(non_sk_utrig_mp_spec,3))');shading flat
% imagesc(lags/Fsd,wfreqs,interp_nonsk_specgram);set(gca,'ydir','normal');
% xlim(xr);
% xlabel('Time (s)');
% ylabel('Frequency (Hz)');
% title('Non-skipped');
% colorbar
% ca = caxis();
% subplot(2,1,2);
% % pcolor(lags/Fsd,wfreqs,squeeze(nanmean(sk_utrig_mp_spec,3))');shading flat
% imagesc(lags/Fsd,wfreqs,interp_sk_specgram);set(gca,'ydir','normal');
% xlim(xr);
% caxis(ca);
% xlabel('Time (s)');
% ylabel('Frequency (Hz)');
% title('Skipped');
% colorbar
% 
% good_hpc_lfp = ~isnan(peak_hpcmua_loc);
% l3mec_pdown_hpclfp = l3mec_pdowns(good_hpc_lfp(l3mec_pdowns));
% used_data = l3mec_pdown_hpclfp;
% 
% non_sk_utrig_hpc_spec = reshape([spec_data(used_data).non_sk_utrig_hpc_spec],[length(lags) length(wfreqs) length(used_data)]);
% sk_utrig_hpc_spec = reshape([spec_data(used_data).sk_utrig_hpc_spec],[length(lags) length(wfreqs) length(used_data)]);
% 
% nonsk_specgram = squeeze(nanmean(non_sk_utrig_hpc_spec,3))';
% interp_nonsk_specgram = interp1(wfreqs,nonsk_specgram,interp_f);
% sk_specgram = squeeze(nanmean(sk_utrig_hpc_spec,3))';
% interp_sk_specgram = interp1(wfreqs,sk_specgram,interp_f);
% 
% h3=figure;
% subplot(2,1,1);
% % pcolor(lags/Fsd,wfreqs,squeeze(nanmean(non_sk_utrig_hpc_spec,3))');shading flat
% imagesc(lags/Fsd,wfreqs,interp_nonsk_specgram);set(gca,'ydir','normal');
% xlim(xr);
% xlabel('Time (s)');
% ylabel('Frequency (Hz)');
% title('Non-skipped');
% ca = caxis();
% colorbar
% subplot(2,1,2);
% % pcolor(lags/Fsd,wfreqs,squeeze(nanmean(sk_utrig_hpc_spec,3))');shading flat
% imagesc(lags/Fsd,wfreqs,interp_sk_specgram);set(gca,'ydir','normal');
% xlim(xr);
% caxis(ca);
% xlabel('Time (s)');
% ylabel('Frequency (Hz)');
% title('Skipped');
% colorbar
% 
% %%
% fig_width = 4.36;
% rel_height = 1.5;
% 
% figufy(h1);
% fname1 = [fig_dir 'utrig_ctx_spec_5s'];
% exportfig(h1,fname1,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(h1);
% 
% figufy(h2);
% fname1 = [fig_dir 'utrig_mp_spec_5s'];
% exportfig(h2,fname1,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(h2);
% 
% figufy(h3);
% fname1 = [fig_dir 'utrig_hpc_spec_5s'];
% exportfig(h3,fname1,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(h3);
% 
% %%
% xr = [-0.5 1.5];
% 
% used_data = l3mec_pups;
% 
% non_sk_dtrig_lf8_spec = reshape([spec_data(used_data).non_sk_dtrig_lfp_spec],[length(lags) length(wfreqs) length(used_data)]);
% sk_dtrig_lf8_spec = reshape([spec_data(used_data).sk_dtrig_lfp_spec],[length(lags) length(wfreqs) length(used_data)]);
% 
% nonsk_specgram = squeeze(nanmean(non_sk_dtrig_lf8_spec,3))';
% interp_nonsk_specgram = interp1(wfreqs,nonsk_specgram,interp_f);
% sk_specgram = squeeze(nanmean(sk_dtrig_lf8_spec,3))';
% interp_sk_specgram = interp1(wfreqs,sk_specgram,interp_f);
% 
% h1=figure;
% subplot(2,1,1);
% % pcolor(lags/Fsd,wfreqs,squeeze(nanmean(non_sk_dtrig_lf8_spec,3))');shading flat
% imagesc(lags/Fsd,wfreqs,interp_nonsk_specgram);set(gca,'ydir','normal');
% xlim(xr);
% xlabel('Time (s)');
% ylabel('Frequency (Hz)');
% title('Non-skipped');
% ca = caxis();
% colorbar
% subplot(2,1,2);
% % pcolor(lags/Fsd,wfreqs,squeeze(nanmean(sk_dtrig_lf8_spec,3))');shading flat
% imagesc(lags/Fsd,wfreqs,interp_sk_specgram);set(gca,'ydir','normal');
% caxis(ca);
% xlim(xr);
% xlabel('Time (s)');
% ylabel('Frequency (Hz)');
% title('Skipped');
% colorbar
% 
% used_data = l3mec_pups;
% 
% non_sk_dtrig_mp_spec = reshape([spec_data(used_data).non_sk_dtrig_mp_spec],[length(lags) length(wfreqs) length(used_data)]);
% sk_dtrig_mp_spec = reshape([spec_data(used_data).sk_dtrig_mp_spec],[length(lags) length(wfreqs) length(used_data)]);
% 
% nonsk_specgram = squeeze(nanmean(non_sk_dtrig_mp_spec,3))';
% interp_nonsk_specgram = interp1(wfreqs,nonsk_specgram,interp_f);
% sk_specgram = squeeze(nanmean(sk_dtrig_mp_spec,3))';
% interp_sk_specgram = interp1(wfreqs,sk_specgram,interp_f);
% 
% h2=figure;
% subplot(2,1,1);
% % pcolor(lags/Fsd,wfreqs,squeeze(nanmean(non_sk_dtrig_mp_spec,3))');shading flat
% imagesc(lags/Fsd,wfreqs,interp_nonsk_specgram);set(gca,'ydir','normal');
% xlim(xr);
% xlabel('Time (s)');
% ylabel('Frequency (Hz)');
% title('Non-skipped');
% ca = caxis();
% colorbar
% subplot(2,1,2);
% % pcolor(lags/Fsd,wfreqs,squeeze(nanmean(sk_dtrig_mp_spec,3))');shading flat
% imagesc(lags/Fsd,wfreqs,interp_sk_specgram);set(gca,'ydir','normal');
% xlim(xr);
% xlabel('Time (s)');
% ylabel('Frequency (Hz)');
% title('Skipped');
% caxis(ca);
% colorbar
% 
% good_hpc_lfp = ~isnan(peak_hpcmua_loc);
% l3mec_pup_hpclfp = l3mec_pups(good_hpc_lfp(l3mec_pups));
% used_data = l3mec_pup_hpclfp;
% 
% non_sk_dtrig_hpc_spec = reshape([spec_data(used_data).non_sk_dtrig_hpc_spec],[length(lags) length(wfreqs) length(used_data)]);
% sk_dtrig_hpc_spec = reshape([spec_data(used_data).sk_dtrig_hpc_spec],[length(lags) length(wfreqs) length(used_data)]);
% 
% nonsk_specgram = squeeze(nanmean(non_sk_dtrig_hpc_spec,3))';
% interp_nonsk_specgram = interp1(wfreqs,nonsk_specgram,interp_f);
% sk_specgram = squeeze(nanmean(sk_dtrig_hpc_spec,3))';
% interp_sk_specgram = interp1(wfreqs,sk_specgram,interp_f);
% 
% h3=figure;
% subplot(2,1,1);
% % pcolor(lags/Fsd,wfreqs,squeeze(nanmean(non_sk_dtrig_hpc_spec,3))');shading flat
% imagesc(lags/Fsd,wfreqs,interp_nonsk_specgram);set(gca,'ydir','normal');
% xlim(xr);
% xlabel('Time (s)');
% ylabel('Frequency (Hz)');
% title('Non-skipped');
% ca = caxis();
% colorbar
% subplot(2,1,2);
% % pcolor(lags/Fsd,wfreqs,squeeze(nanmean(sk_dtrig_hpc_spec,3))');shading flat
% imagesc(lags/Fsd,wfreqs,interp_sk_specgram);set(gca,'ydir','normal');
% xlim(xr);
% xlabel('Time (s)');
% ylabel('Frequency (Hz)');
% title('Skipped');
% caxis(ca);
% colorbar
% 
% %%
% fig_width = 4.36;
% rel_height = 1.5;
% 
% figufy(h1);
% fname1 = [fig_dir 'dtrig_ctx_spec_5s'];
% exportfig(h1,fname1,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(h1);
% 
% figufy(h2);
% fname1 = [fig_dir 'dtrig_mp_spec_5s'];
% exportfig(h2,fname1,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(h2);
% 
% figufy(h3);
% fname1 = [fig_dir 'dtrig_hpc_spec_5s'];
% exportfig(h3,fname1,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(h3);
% 
% close all
% 
% %%
% close all
% h1=figure;
% subplot(2,1,1);
% pcolor(lags/Fsd,wfreqs,squeeze(nanmean(sk_utrig_mp_spec,3))');shading flat
% xlim(xr);
% xlabel('Time (s)');
% ylabel('Frequency (Hz)');
% title('Skipped up-triggered MP');
% caxis([-1 -0.25])
% colorbar
% subplot(2,1,2);
% pcolor(lags/Fsd,wfreqs,squeeze(nanmean(sk_dtrig_mp_spec,3))');shading flat
% xlim(xr);
% xlabel('Time (s)');
% ylabel('Frequency (Hz)');
% title('Skipped');
% caxis([0 0.75]);
% colorbar
% title('Skipped down-triggered MP');
% 
% fig_width = 4.36;
% rel_height = 1.5;
% 
% figufy(h1);
% fname1 = [fig_dir 'sk_trig_mp_spec'];
% exportfig(h1,fname1,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(h1);
% 
% %%
% l3mec_non_sk_utrig_mp_spec = reshape([spec_data(l3mec).non_sk_utrig_mp_spec],[length(lags) length(wfreqs) length(l3mec)]);
% l3mec_avg_mp_spec = reshape([spec_data(l3mec).avg_mp_spectrum],[length(wfreqs) length(l3mec)]);
% l3mec_std_mp_spec = reshape([spec_data(l3mec).std_mp_spectrum],[length(wfreqs) length(l3mec)]);
% 
% used_l3mec = find(max(l3mec_avg_mp_spec) < 1);
% 
% % l3mec_non_sk_utrig_mp_spec = bsxfun(@times,l3mec_non_sk_utrig_mp_spec,reshape(l3mec_std_mp_spec,[1 length(wfreqs) length(l3mec)]));
% % l3mec_non_sk_utrig_mp_spec = bsxfun(@plus,l3mec_non_sk_utrig_mp_spec,reshape(l3mec_avg_mp_spec,[1 length(wfreqs) length(l3mec)]));
% 
% l3lec_non_sk_utrig_mp_spec = reshape([spec_data(l3lec).non_sk_utrig_mp_spec],[length(lags) length(wfreqs) length(l3lec)]);
% l3lec_avg_mp_spec = reshape([spec_data(l3lec).avg_mp_spectrum],[length(wfreqs) length(l3lec)]);
% l3lec_std_mp_spec = reshape([spec_data(l3lec).std_mp_spectrum],[length(wfreqs) length(l3lec)]);
% 
% used_l3lec = find(max(l3lec_avg_mp_spec) < 1);
% 
% 
% % l3lec_non_sk_utrig_mp_spec = bsxfun(@times,l3lec_non_sk_utrig_mp_spec,reshape(l3lec_std_mp_spec,[1 length(wfreqs) length(l3lec)]));
% % l3lec_non_sk_utrig_mp_spec = bsxfun(@plus,l3lec_non_sk_utrig_mp_spec,reshape(l3lec_avg_mp_spec,[1 length(wfreqs) length(l3lec)]));
% 
% % figure
% % subplot(2,1,1)
% % pcolor(lags/Fsd,wfreqs,squeeze(nanmean(l3mec_non_sk_utrig_mp_spec(:,:,used_l3mec),3))');shading flat
% % xlim([-1 1.5]);
% % subplot(2,1,2)
% % pcolor(lags/Fsd,wfreqs,squeeze(nanmean(l3lec_non_sk_utrig_mp_spec(:,:,used_l3lec),3))');shading flat
% % xlim([-1 1.5]);
% 
% pow_diff = squeeze(nanmean(l3mec_non_sk_utrig_mp_spec(:,:,used_l3mec),3))' - squeeze(nanmean(l3lec_non_sk_utrig_mp_spec(:,:,used_l3lec),3))';
% 
% figure
% pcolor(lags/Fsd,wfreqs,pow_diff);shading flat
% xlim([-1 1.5]);
% 
% %%
% l3mec_non_sk_dtrig_mp_spec = reshape([spec_data(l3mec).non_sk_dtrig_mp_spec],[length(lags) length(wfreqs) length(l3mec)]);
% l3mec_avg_mp_spec = reshape([spec_data(l3mec).avg_mp_spectrum],[length(wfreqs) length(l3mec)]);
% l3mec_std_mp_spec = reshape([spec_data(l3mec).std_mp_spectrum],[length(wfreqs) length(l3mec)]);
% 
% used_l3mec = find(max(l3mec_avg_mp_spec) < 1);
% 
% % l3mec_non_sk_utrig_mp_spec = bsxfun(@times,l3mec_non_sk_utrig_mp_spec,reshape(l3mec_std_mp_spec,[1 length(wfreqs) length(l3mec)]));
% % l3mec_non_sk_utrig_mp_spec = bsxfun(@plus,l3mec_non_sk_utrig_mp_spec,reshape(l3mec_avg_mp_spec,[1 length(wfreqs) length(l3mec)]));
% 
% l3lec_non_sk_dtrig_mp_spec = reshape([spec_data(l3lec).non_sk_dtrig_mp_spec],[length(lags) length(wfreqs) length(l3lec)]);
% l3lec_avg_mp_spec = reshape([spec_data(l3lec).avg_mp_spectrum],[length(wfreqs) length(l3lec)]);
% l3lec_std_mp_spec = reshape([spec_data(l3lec).std_mp_spectrum],[length(wfreqs) length(l3lec)]);
% 
% used_l3lec = find(max(l3lec_avg_mp_spec) < 1);
% 
% 
% % l3lec_non_sk_utrig_mp_spec = bsxfun(@times,l3lec_non_sk_utrig_mp_spec,reshape(l3lec_std_mp_spec,[1 length(wfreqs) length(l3lec)]));
% % l3lec_non_sk_utrig_mp_spec = bsxfun(@plus,l3lec_non_sk_utrig_mp_spec,reshape(l3lec_avg_mp_spec,[1 length(wfreqs) length(l3lec)]));
% 
% pow_diff = squeeze(nanmean(l3mec_non_sk_dtrig_mp_spec(:,:,used_l3mec),3))' - squeeze(nanmean(l3lec_non_sk_dtrig_mp_spec(:,:,used_l3lec),3))';
% 
% figure
% pcolor(lags/Fsd,wfreqs,pow_diff);shading flat
% xlim([-1 1.5]);
