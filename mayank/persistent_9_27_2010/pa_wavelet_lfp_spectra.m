clear all

load F:\WC_Germany\overall_EC\overall_EC_dir.mat
addpath('F:\WC_Germany\parietal_cortical_2010\')
addpath('F:\WC_Germany\hsmm_state_detection\')
addpath('F:\WC_Germany\persistent_2010\')
addpath('F:\Code\smoothing\software\')
addpath('F:\Code\general\')
addpath('F:\Code\wavelet_tools')

used_data = [l3mec_p l3lec_p];
sess_data = sess_data(used_data);

Fs = 2016;
dsf = 8;
Fsd = 2016/dsf;
niqf = 2016/2;
[b,a] = butter(2,[0.05/niqf 120/niqf]);
[b_lf,a_lf] = butter(2,[0.05/niqf 40/niqf]);
[b_hf,a_hf] = butter(2,[15/niqf 80/niqf]);
hf_smooth = round(Fsd*0.05);

up_range = [0.3 7];
down_range = [0.3 7];
numBins = 25;
lin_grid = linspace(up_range(1),up_range(2),numBins+1);
log_grid = logspace(log10(up_range(1)),log10(up_range(2)),numBins+1);

edge_dur = 0.1;
min_state_dur = 0.5;

% Fs = 2016;
% dsf = 16;
% Fsd = 2016/dsf;
% niqf = 2016/2;
% [b,a] = butter(2,[0.05/niqf 50/niqf]);
% 

%% CWT params

min_freq = 0.4; max_freq = 100; delta_j = 0.15;
k0 = 6; %wavenumber for morlet wavelet
fourier_factor = (4*pi)/(k0 + sqrt(2 + k0^2));
min_scale = 1/max_freq/fourier_factor;
max_scale = 1/min_freq/fourier_factor;
n_scales = round(1/delta_j*log2(max_scale/min_scale));

%%
nlx_amp_range = linspace(-4,4,100);

load ./pa_corresponding_lfp_state_data_rtest

for d = 1:length(sess_data)
    cdir = sess_data(d).directory;
    cdir(1) = 'F';
    disp(sprintf('session %d',d))
    cd(cdir);
    s_name = strcat(sess_data(d).region,'_l',sess_data(d).layer(1),'_',sess_data(d).name);
    
    load ./used_data lf8 wcv_minus_spike lf3
    lf8_d = zscore(downsample(filtfilt(b,a,lf8),dsf));
    wcv_d = zscore(downsample(filtfilt(b,a,wcv_minus_spike),dsf));
    lf3_d = zscore(downsample(filtfilt(b,a,lf3),dsf));
        
    lf8_lf = zscore(downsample(filtfilt(b_lf,a_lf,lf8),dsf));
%     wcv_lf = zscore(downsample(filtfilt(b_lf,a_lf,wcv_minus_spike),dsf));
%     lf3_lf = zscore(downsample(filtfilt(b_lf,a_lf,lf3),dsf));
    
    lf8_hf = zscore(downsample(filtfilt(b_hf,a_hf,lf8),dsf));
    wcv_hf = zscore(downsample(filtfilt(b_hf,a_hf,wcv_minus_spike),dsf));
    lf3_hf = zscore(downsample(filtfilt(b_hf,a_hf,lf3),dsf));

    lf8_hf = sqrt(jmm_smooth_1d_cor(lf8_hf.^2,hf_smooth));
    wcv_hf = sqrt(jmm_smooth_1d_cor(wcv_hf.^2,hf_smooth));
    lf3_hf = sqrt(jmm_smooth_1d_cor(lf3_hf.^2,hf_smooth));

    time = (1:length(wcv_d))/Fsd;
    time = time(:);
        
    load ./pa_hsmm_state_seq
    mp_state_seq = hsmm_bbstate_seq;
    [new_seg_inds] = resample_uds_seg_inds_v2(hmm.UDS_segs,hmm.Fs,252,mp_state_seq);
    [up_trans_inds,down_trans_inds] = compute_state_transitions_seg(new_seg_inds,mp_state_seq);
    c_dsf = 252/Fsd;
    up_trans_inds = round(up_trans_inds/c_dsf);
    down_trans_inds = round(down_trans_inds/c_dsf);
    
    load ./pa_hsmm_state_seq8.mat
    lfp_state_seq = hsmm_bbstate_seq8;
    [new_seg_inds8] = resample_uds_seg_inds_v2(hmm8.UDS_segs,hmm8.Fs,252,lfp_state_seq);
    [up_trans_inds8,down_trans_inds8] = compute_state_transitions_seg(new_seg_inds8,lfp_state_seq);
    up_trans_inds8 = round(up_trans_inds8/c_dsf);
    down_trans_inds8 = round(down_trans_inds8/c_dsf);
    
    [lfp_state_durations] = compute_state_durations_seg(lfp_state_seq,Fsd);
    [mp_state_durations] = compute_state_durations_seg(mp_state_seq,Fsd);
    
    d_new_seg_inds = ceil(new_seg_inds/c_dsf);
    uds_inds = zeros(size(time));
    for i = 1:size(d_new_seg_inds,1)
        uds_inds(d_new_seg_inds(i,1):d_new_seg_inds(i,2)) = 1;
    end

    mp_state_vec = zeros(size(time));
    for i = 1:length(up_trans_inds)
        mp_state_vec(up_trans_inds(i):down_trans_inds(i))=1;
    end
    lf8_state_vec = zeros(size(time));
    for i = 1:length(up_trans_inds8)
        lf8_state_vec(up_trans_inds8(i):down_trans_inds8(i))=1;
    end
    mp_state_vec(uds_inds==0) = nan;
    lf8_state_vec(uds_inds==0) = nan;
    %     mp_state_vec = logical(mp_state_vec);
    %     lf8_state_vec = logical(lf8_state_vec);
    
%     lf8_long_vec = nan(size(time));
    lf8_skipped_vec = nan(size(time));
    lf8_upstate_pow = nan(size(down_trans_inds8));
    lf8_downstate_pow = nan(size(down_trans_inds8));
    lf8_upstate_pow_ne = nan(size(down_trans_inds8));
    lf8_downstate_pow_ne = nan(size(down_trans_inds8));
    for i = 1:length(down_trans_inds8)-1
        curdownset = down_trans_inds8(i):up_trans_inds8(i+1);
        curupset = up_trans_inds8(i):down_trans_inds8(i);
        lf8_upstate_pow(i) = mean(lf8_hf(curupset));
        lf8_downstate_pow(i) = mean(lf8_hf(curdownset));
        if length(curdownset)/Fsd > min_state_dur-2*edge_dur
            %get rid of transition regions
            curdownset(1:round(Fsd*edge_dur)) = [];
            curdownset(end-round(edge_dur)+1:end) = [];
            lf8_downstate_pow_ne(i) = mean(lf8_hf(curdownset));
        end
        if length(curupset)/Fsd > min_state_dur-2*edge_dur
            curupset(1:round(Fsd*edge_dur)) = [];
            curupset(end-round(Fsd*edge_dur)+1:end) = [];
            lf8_upstate_pow_ne(i) = mean(lf8_hf(curupset));
        end
        if length(curdownset)/Fsd > 0.5-2*edge_dur %only use down states of sufficient duration
            if all(mp_state_vec(curdownset))
                lf8_skipped_vec(curdownset) = 1;
            else
                lf8_skipped_vec(curdownset) = 0;
            end
%             lf8_long_vec(curdownset) = 0;
        end
%         if length(curdownset)/Fsd > 0.8 
%             lf8_long_vec(curdownset) = 1;
%         end
    end
    lf8_skipped_vec(uds_inds==0) = nan;
%     lf8_long_vec(uds_inds==0) = nan;
    
%     lf8_trig_vec = zeros(size(time));
%     trig_up8_ids = corresp_lfp_upid{d}(rpersd_across_ids{d});
%     for i = 1:length(trig_up8_ids)
%         curupset = up_trans_inds8(trig_up8_ids(i)):down_trans_inds8(trig_up8_ids(i));
%         lf8_trig_vec(curupset) = 1;
%     end
%     lf8_trig_vec(uds_inds==0) = nan;
        
   
    mp_upstate_pow = nan(size(down_trans_inds));
    mp_downstate_pow = nan(size(down_trans_inds));
    mp_upstate_pow_ne = nan(size(down_trans_inds));
    mp_downstate_pow_ne = nan(size(down_trans_inds));
    for i = 1:length(down_trans_inds)-1
        curdownset = down_trans_inds(i):up_trans_inds(i+1);
        curupset = up_trans_inds(i):down_trans_inds(i);
        mp_upstate_pow(i) = mean(wcv_hf(curupset));
        mp_downstate_pow(i) = mean(wcv_hf(curdownset));
        if length(curdownset)/Fsd > min_state_dur-2*edge_dur
            curdownset(1:round(edge_dur*Fsd)) = [];
            curdownset(end-round(edge_dur*Fsd)+1:end) = [];
            mp_downstate_pow_ne(i) = mean(wcv_hf(curdownset));
        end
        if length(curupset)/Fsd > min_state_dur-2*edge_dur
            curupset(1:round(edge_dur*Fsd)) = [];
            curupset(end-round(edge_dur*Fsd)+1:end) = [];
            mp_upstate_pow_ne(i) = mean(wcv_hf(curupset));
        end
    end
    mp_updur_lfpc_delay{d}(length(down_trans_inds):end) = nan;
    
    lfp_dskipped_downs{d}(lfp_state_durations{1}(lfp_dskipped_downs{d}) < min_state_dur) = [];
    lfp_ndskipped_downs{d}(lfp_state_durations{1}(lfp_ndskipped_downs{d}) < min_state_dur) = [];
    lfp_skipped_downs{d}(lfp_state_durations{1}(lfp_skipped_downs{d}) < min_state_dur) = [];
    lfp_nskipped_downs{d}(lfp_state_durations{1}(lfp_nskipped_downs{d}) < min_state_dur) = [];
        
    mp_t2pers_pow(d) = nanmean(mp_upstate_pow_ne(rpers_across_ids{d}));
    mp_nt2pers_pow(d) = nanmean(mp_upstate_pow_ne(mp_updur_lfpc_delay{d} < 1));

%     lfp_dskippeddown_pow(d) = nanmean(lf8_downstate_pow_ne(lfp_dskipped_downs{d}));
%     lfp_ndskippeddown_pow(d) = nanmean(lf8_downstate_pow_ne(lfp_ndskipped_downs{d}));
    lfp_dskippeddown_pow(d) = nanmean(lf8_downstate_pow(lfp_dskipped_downs{d}));
    lfp_ndskippeddown_pow(d) = nanmean(lf8_downstate_pow(lfp_ndskipped_downs{d}));
    lfp_skippeddown_pow(d) = nanmean(lf8_downstate_pow(lfp_skipped_downs{d}));
    lfp_nskippeddown_pow(d) = nanmean(lf8_downstate_pow(lfp_nskipped_downs{d}));
    
    lfp_trigup_pow(d) = nanmean(lf8_upstate_pow_ne(corresp_lfp_upid{d}(rpers_across_ids{d})));
    lfp_ntrigup_pow(d) = nanmean(lf8_upstate_pow_ne(corresp_lfp_upid{d}(mp_updur_lfpc_delay{d} < 1)));
    
    
    lfp_ndskipped_downs{d}(lfp_ndskipped_downs{d}==length(down_trans_inds8)) = [];
    lfp_nskipped_downs{d}(lfp_nskipped_downs{d}==length(down_trans_inds8)) = [];
    
    lfp_nskip_down_vec = zeros(size(time));
    for i = 1:length(lfp_nskipped_downs{d})
        curdownset = down_trans_inds8(lfp_nskipped_downs{d}(i)):up_trans_inds8(lfp_nskipped_downs{d}(i)+1);
        curdownset(1:round(Fsd*edge_dur)) = [];
        curdownset(end-round(edge_dur*Fsd)+1:end) = [];
        lfp_nskip_down_vec(curdownset) = 1;
%         curdownset = curdownset + round(lfp_up_lag{d}(i)*Fsd);
%         lfp_dnskip_down_vec(curdownset) = 1;
    end
    lfp_t1_pow(d) = mean(lf8_hf(lfp_nskip_down_vec==1 & mp_state_vec==1));
    lfp_nt1_pow(d) = mean(lf8_hf(lfp_nskip_down_vec==1 & mp_state_vec==0));
    
    mp_upstate_durpow(d,:) = nan(size(log_grid));
    mp_downstate_durpow(d,:) = nan(size(log_grid));
    mp_upstate_durpow_ne(d,:) = nan(size(log_grid));
    mp_downstate_durpow_ne(d,:) = nan(size(log_grid));
    lf8_upstate_durpow(d,:) = nan(size(log_grid));
    lf8_downstate_durpow(d,:) = nan(size(log_grid));
    lf8_upstate_durpow_ne(d,:) = nan(size(log_grid));
    lf8_downstate_durpow_ne(d,:) = nan(size(log_grid));
   for j = 1:length(log_grid)-1
        curset = find(mp_state_durations{2} >= log_grid(j) & mp_state_durations{2} < log_grid(j+1));
        if ~isempty(curset)
            mp_upstate_durpow(d,j) = nanmean(mp_upstate_pow(curset));
            mp_upstate_durpow_ne(d,j) = nanmean(mp_upstate_pow_ne(curset));
        end
        curset = find(mp_state_durations{1} >= log_grid(j) & mp_state_durations{1} < log_grid(j+1));
        if ~isempty(curset)
            mp_downstate_durpow(d,j) = nanmean(mp_downstate_pow(curset));
            mp_downstate_durpow_ne(d,j) = nanmean(mp_downstate_pow_ne(curset));
        end
        curset = find(lfp_state_durations{2} >= log_grid(j) & lfp_state_durations{2} < log_grid(j+1));
        if ~isempty(curset)
            lf8_upstate_durpow(d,j) = nanmean(lf8_upstate_pow(curset));
            lf8_upstate_durpow_ne(d,j) = nanmean(lf8_upstate_pow_ne(curset));
        end
        curset = find(lfp_state_durations{1} >= log_grid(j) & lfp_state_durations{1} < log_grid(j+1));
        if ~isempty(curset)
            lf8_downstate_durpow(d,j) = nanmean(lf8_downstate_pow(curset));
            lf8_downstate_durpow_ne(d,j) = nanmean(lf8_downstate_pow_ne(curset));
        end        
    end
    
%     [wcv_scalogram,periods,scales] = wavelet(wcv_d,1/Fsd,1,delta_j,min_scale,n_scales);
%     cwt_freqs = 1./periods;    
%     [lf8_scalogram,periods,scales] = wavelet(lf8_d,1/Fsd,1,delta_j,min_scale,n_scales);
%     [lf3_scalogram,periods,scales] = wavelet(lf3_d,1/Fsd,1,delta_j,min_scale,n_scales);
%     
%     inv_scales = 1./scales';
%     n = size(wcv_scalogram,2);
% %     sm_wcv_scalogram = smoothwavelet(inv_scales(:,ones(1,n)).*abs(wcv_scalogram).^2,1,delta_j,scales);
% %     sm_lf8_scalogram = smoothwavelet(inv_scales(:,ones(1,n)).*abs(lf8_scalogram).^2,1,delta_j,scales);
% %     sm_lf3_scalogram = smoothwavelet(inv_scales(:,ones(1,n)).*abs(lf3_scalogram).^2,1,delta_j,scales);
%     
% %     xw8_scalogram=wcv_scalogram.*conj(lf8_scalogram);
% %     xw3_scalogram=wcv_scalogram.*conj(lf3_scalogram);
%     
% %     s_xw8_scalogram = smoothwavelet(inv_scales(:,ones(1,n)).*xw8_scalogram,1,delta_j,scales);
% %     s_xw3_scalogram = smoothwavelet(inv_scales(:,ones(1,n)).*xw3_scalogram,1,delta_j,scales);
%     
% %     cw8_coherence = abs(s_xw8_scalogram).^2./(sm_wcv_scalogram.*sm_lf8_scalogram);
% %     cw3_coherence = abs(s_xw3_scalogram).^2./(sm_wcv_scalogram.*sm_lf3_scalogram);
%     
%     %% to log transform
%     wcv_scalogram = 10*log10(abs(wcv_scalogram).^2);
%     lf8_scalogram = 10*log10(abs(lf8_scalogram).^2);
%     lf3_scalogram = 10*log10(abs(lf3_scalogram).^2);
%     
%     %% to z-score amplitudes separately within each state
% %     wcv_scalogram = abs(wcv_scalogram);
% %     lf8_scalogram = abs(lf8_scalogram);
% %     lf3_scalogram = abs(lf3_scalogram);
% 
%     wcv_scalogram_upz = wcv_scalogram - repmat(mean(wcv_scalogram(:,mp_state_vec==1),2),1,length(time));
%     wcv_scalogram_upz = wcv_scalogram_upz./repmat(std(wcv_scalogram(:,mp_state_vec==1),[],2),1,length(time));
%     wcv_scalogram_downz = wcv_scalogram - repmat(mean(wcv_scalogram(:,mp_state_vec==0),2),1,length(time));
%     wcv_scalogram_downz = wcv_scalogram_downz./repmat(std(wcv_scalogram(:,mp_state_vec==0),[],2),1,length(time));
% %     
%     lf8_scalogram_upz = lf8_scalogram - repmat(mean(lf8_scalogram(:,lf8_state_vec==1),2),1,length(time));
%     lf8_scalogram_upz = lf8_scalogram_upz./repmat(std(lf8_scalogram(:,lf8_state_vec==1),[],2),1,length(time));
%     lf8_scalogram_downz = lf8_scalogram - repmat(mean(lf8_scalogram(:,lf8_state_vec==0),2),1,length(time));
%     lf8_scalogram_downz = lf8_scalogram_downz./repmat(std(lf8_scalogram(:,lf8_state_vec==0),[],2),1,length(time));
% 
%     wcv_scalogram = zscore(wcv_scalogram')';
%     lf8_scalogram = zscore(lf8_scalogram')';
%     lf3_scalogram = zscore(lf3_scalogram')';
%     
% %     mp_up_spec(d,:) = nanmean(wcv_scalogram(:,mp_state_vec==1),2);
% %     mp_down_spec(d,:) = nanmean(wcv_scalogram(:,mp_state_vec==0),2);
%     mp_up_up8_spec(d,:) = nanmean(wcv_scalogram_upz(:,mp_state_vec==1 & lf8_state_vec==1),2);
%     mp_up_down8_spec(d,:) = nanmean(wcv_scalogram_upz(:,mp_state_vec==1 & lf8_state_vec==0),2);
%     mp_down_up8_spec(d,:) = nanmean(wcv_scalogram_downz(:,mp_state_vec==0 & lf8_state_vec==1),2);
%     mp_down_down8_spec(d,:) = nanmean(wcv_scalogram_downz(:,mp_state_vec==0 & lf8_state_vec==0),2);
%     
% %     lf8_up8_spec(d,:) = nanmean(lf8_scalogram(:,lf8_state_vec==1),2);
% %     lf8_down8_spec(d,:) = nanmean(lf8_scalogram(:,lf8_state_vec==0),2);
%     lf8_up_up8_spec(d,:) = nanmean(lf8_scalogram_upz(:,mp_state_vec==1 & lf8_state_vec==1),2);
%     lf8_up_down8_spec(d,:) = nanmean(lf8_scalogram_downz(:,mp_state_vec==1 & lf8_state_vec==0),2);
%     lf8_down_up8_spec(d,:) = nanmean(lf8_scalogram_upz(:,mp_state_vec==0 & lf8_state_vec==1),2);
%     lf8_down_down8_spec(d,:) = nanmean(lf8_scalogram_downz(:,mp_state_vec==0 & lf8_state_vec==0),2);
%     
%     if nansum(lf8_skipped_vec) >= round(Fsd)
%         lf8_down_skipped_spec(d,:) = nanmean(lf8_scalogram_downz(:,lf8_skipped_vec==1),2);
%         lf8_down_nonskipped_spec(d,:) = nanmean(lf8_scalogram_downz(:,lf8_skipped_vec==0),2);
%         lf8_down_skipped_amp(d,:) = ksdensity(lf8_lf(lf8_skipped_vec==1),nlx_amp_range);
%         lf8_down_nonskipped_amp(d,:) = ksdensity(lf8_lf(lf8_skipped_vec==0),nlx_amp_range);
%     else
%         lf8_down_skipped_spec(d,:) = nan;
%         lf8_down_skipped_amp(d,:) = nan;
%         lf8_down_nonskipped_spec(d,:) = nan;
%         lf8_down_nonskipped_amp(d,:) = nan;
%     end
% 
%     if nansum(lf8_long_vec) >= round(Fsd)
%         lf8_down_long_spec(d,:) = nanmean(lf8_scalogram_downz(:,lf8_long_vec==1),2);
%         lf8_down_nonlong_spec(d,:) = nanmean(lf8_scalogram_downz(:,lf8_long_vec==0),2);
%         lf8_down_long_amp(d,:) = ksdensity(lf8_lf(lf8_long_vec==1),nlx_amp_range);
%         lf8_down_nonlong_amp(d,:) = ksdensity(lf8_lf(lf8_long_vec==0),nlx_amp_range);
%     else
%         lf8_down_long_spec(d,:) = nan;
%         lf8_down_long_amp(d,:) = nan;
%         lf8_down_long_spec(d,:) = nan;
%         lf8_down_long_amp(d,:) = nan;
%     end
%     
% %     lf3_up8_spec(d,:) = nanmean(lf3_scalogram(:,lf8_state_vec==1),2);
% %     lf3_down8_spec(d,:) = nanmean(lf3_scalogram(:,lf8_state_vec==0),2);
%     lf3_up_up8_spec(d,:) = nanmean(lf3_scalogram(:,mp_state_vec==1 & lf8_state_vec==1),2);
%     lf3_up_down8_spec(d,:) = nanmean(lf3_scalogram(:,mp_state_vec==1 & lf8_state_vec==0),2);
%     lf3_down_up8_spec(d,:) = nanmean(lf3_scalogram(:,mp_state_vec==0 & lf8_state_vec==1),2);
%     lf3_down_down8_spec(d,:) = nanmean(lf3_scalogram(:,mp_state_vec==0 & lf8_state_vec==0),2);
%         
% %     cw8_up8_spec(d,:) = nanmean(cw8_coherence(:,lf8_state_vec==1),2);
% %     cw8_down8_spec(d,:) = nanmean(cw8_coherence(:,lf8_state_vec==0),2);
% %     cw8_up_up8_spec(d,:) = nanmean(cw8_coherence(:,mp_state_vec==1 & lf8_state_vec==1),2);
% %     cw8_up_down8_spec(d,:) = nanmean(cw8_coherence(:,mp_state_vec==1 & lf8_state_vec==0),2);
% %     cw8_down_up8_spec(d,:) = nanmean(cw8_coherence(:,mp_state_vec==0 & lf8_state_vec==1),2);
% %     cw8_down_down8_spec(d,:) = nanmean(cw8_coherence(:,mp_state_vec==0 & lf8_state_vec==0),2);
% %        
% %     cw3_up8_spec(d,:) = nanmean(cw3_coherence(:,lf8_state_vec==1),2);
% %     cw3_down8_spec(d,:) = nanmean(cw3_coherence(:,lf8_state_vec==0),2);
% %     cw3_up_up8_spec(d,:) = nanmean(cw3_coherence(:,mp_state_vec==1 & lf8_state_vec==1),2);
% %     cw3_up_down8_spec(d,:) = nanmean(cw3_coherence(:,mp_state_vec==1 & lf8_state_vec==0),2);
% %     cw3_down_up8_spec(d,:) = nanmean(cw3_coherence(:,mp_state_vec==0 & lf8_state_vec==1),2);
% %     cw3_down_down8_spec(d,:) = nanmean(cw3_coherence(:,mp_state_vec==0 & lf8_state_vec==0),2);
       
end

cd F:\WC_Germany\persistent_9_27_2010\
save pa_wvlet_powspectra_robust_z2 mp_* lf8_* skipped* lfp_*

%%
mec_cells = 1:length(l3mec_p);
lec_cells = (length(l3mec_p)+1):(length(l3mec_p)+length(l3lec_p));

%%
figure
% h = errorbar(log_grid,nanmean(lf8_downstate_durpow),nanstd(lf8_downstate_durpow)/sqrt(36));
hold on
% h = errorbar(log_grid,nanmean(lf8_upstate_durpow),nanstd(lf8_upstate_durpow)/sqrt(36),'r');
h = errorbar(log_grid,nanmean(lf8_downstate_durpow_ne),nanstd(lf8_downstate_durpow_ne)/sqrt(36),'k');
h = errorbar(log_grid,nanmean(lf8_upstate_durpow_ne),nanstd(lf8_upstate_durpow_ne)/sqrt(36),'g');
xlim([0 5])
set(gca,'yscale','log')
xlabel('State Duration (s)','Fontsize',14)
ylabel('Average HF Power (AU)','FontSize',14)
legend('Down state','Up state')

figure
% h = errorbar(log_grid,nanmean(mp_downstate_durpow),nanstd(mp_downstate_durpow)/sqrt(36));
hold on
% h = errorbar(log_grid,nanmean(mp_upstate_durpow),nanstd(mp_upstate_durpow)/sqrt(36),'r');
h = errorbar(log_grid,nanmean(mp_downstate_durpow_ne(mec_cells,:)),nanstd(mp_downstate_durpow_ne(mec_cells,:))/sqrt(22),'k');
h = errorbar(log_grid,nanmean(mp_upstate_durpow_ne(mec_cells,:)),nanstd(mp_upstate_durpow_ne(mec_cells,:))/sqrt(22),'g');
% h = errorbar(log_grid,nanmean(mp_downstate_durpow_ne(lec_cells,:)),nanstd(mp_downstate_durpow_ne(lec_cells,:))/sqrt(14),'b');
% h = errorbar(log_grid,nanmean(mp_upstate_durpow_ne(lec_cells,:)),nanstd(mp_upstate_durpow_ne(lec_cells,:))/sqrt(14),'r');
xlim([0 5])

%%
pred_skippeddown_durpow = nansum(lf8_downstate_durpow_ne.*skipped_down_durs,2);
pred_nskippeddown_durpow = nansum(lf8_downstate_durpow_ne.*nonskipped_down_durs,2);

figure
plot(lfp_dskippeddown_pow(mec_cells),lfp_ndskippeddown_pow(mec_cells),'o')
hold on
plot(lfp_skippeddown_pow(mec_cells),lfp_nskippeddown_pow(mec_cells),'ro')
% plot(lfp_dskippeddown_pow(mec_cells),pred_skippeddown_durpow(mec_cells),'ro')
line([0 1],[0 1],'color','k')
xlabel('Skipped Down state HF Power (AU)','FOntsize',14)
ylabel('Non-skipped Down state HF Power (AU)','FOntsize',14)
xlim([0.25 1])
ylim([0.25 1])

figure
plot(lfp_trigup_pow(mec_cells),lfp_ntrigup_pow(mec_cells),'o')
hold on
% plot(lfp_dskippeddown_pow(mec_cells),pred_skippeddown_durpow(mec_cells),'ro')
line([1 2.2],[1 2.2],'color','k')
xlabel('Skipped Down state HF Power (AU)','FOntsize',14)
ylabel('Non-skipped Down state HF Power (AU)','FOntsize',14)
xlim([1 2.2])
ylim([1 2.2])

figure
plot(mp_t2pers_pow(mec_cells),mp_nt2pers_pow(mec_cells),'o')
hold on
% plot(lfp_dskippeddown_pow(mec_cells),pred_skippeddown_durpow(mec_cells),'ro')
line([0.8 2.],[0.8 2.],'color','k')
xlabel('Persistent UP state HF Power (AU)','FOntsize',14)
ylabel('Non-persistent Up state HF Power (AU)','FOntsize',14)
xlim([0.8 2.])
ylim([0.8 2.])

% figure
% plot(lfp_dskippeddown_pow(lec_cells),lfp_ndskippeddown_pow(lec_cells),'o')
% hold on
% plot(lfp_dskippeddown_pow(lec_cells),pred_skippeddown_durpow(lec_cells),'ro')
% line([0 1],[0 1],'color','k')
%%
% figure
% h = errorbar(cwt_freqs,nanmean(mp_up_spec(mec_cells,:)),nanstd(mp_up_spec(mec_cells,:))/sqrt(length(mec_cells)),'b');
% hold on
% h = errorbar(cwt_freqs,nanmean(mp_down_spec(mec_cells,:)),nanstd(mp_up_spec(mec_cells,:))/sqrt(length(mec_cells)),'r');
% h = errorbar(cwt_freqs,nanmean(mp_up_up8_spec(mec_cells,:)),nanstd(mp_up_up8_spec(mec_cells,:))/sqrt(length(mec_cells)),'k');
% h = errorbar(cwt_freqs,nanmean(mp_up_down8_spec(mec_cells,:)),nanstd(mp_up_down8_spec(mec_cells,:))/sqrt(length(mec_cells)),'g');
% h = errorbar(cwt_freqs,nanmean(mp_down_up8_spec(mec_cells,:)),nanstd(mp_down_up8_spec(mec_cells,:))/sqrt(length(mec_cells)),'k');
% h = errorbar(cwt_freqs,nanmean(mp_down_down8_spec(mec_cells,:)),nanstd(mp_down_down8_spec(mec_cells,:))/sqrt(length(mec_cells)),'g');
% set(gca,'xscale','log')
% xlim([2 100])
% set(gca,'FontName','Arial','Fontsize',12)
% xlabel('Frequency (Hz)','Fontsize',12,'fontname','Arial')
% ylabel('Relative Amplitude (AU)','Fontsize',12,'fontname','Arial')
% 
% figure
% h = errorbar(cwt_freqs,nanmean(mp_up_spec(lec_cells,:)),nanstd(mp_up_spec(lec_cells,:))/sqrt(length(lec_cells)),'b');
% hold on
% h = errorbar(cwt_freqs,nanmean(mp_down_spec(lec_cells,:)),nanstd(mp_up_spec(lec_cells,:))/sqrt(length(lec_cells)),'r');
% h = errorbar(cwt_freqs,nanmean(mp_up_up8_spec(lec_cells,:)),nanstd(mp_up_up8_spec(lec_cells,:))/sqrt(length(lec_cells)),'k');
% h = errorbar(cwt_freqs,nanmean(mp_up_down8_spec(lec_cells,:)),nanstd(mp_up_down8_spec(lec_cells,:))/sqrt(length(lec_cells)),'g');
% h = errorbar(cwt_freqs,nanmean(mp_down_up8_spec(lec_cells,:)),nanstd(mp_down_up8_spec(lec_cells,:))/sqrt(length(lec_cells)),'k');
% h = errorbar(cwt_freqs,nanmean(mp_down_down8_spec(lec_cells,:)),nanstd(mp_down_down8_spec(lec_cells,:))/sqrt(length(lec_cells)),'g');
% set(gca,'xscale','log')
% xlim([2 100])
% set(gca,'FontName','Arial','Fontsize',12)
% xlabel('Frequency (Hz)','Fontsize',12,'fontname','Arial')
% ylabel('Relative Amplitude (AU)','Fontsize',12,'fontname','Arial')

mp_up_diff8 = mp_up_up8_spec - mp_up_down8_spec;
mp_down_diff8 = mp_down_up8_spec - mp_down_down8_spec;

figure
h = errorbar(cwt_freqs,nanmean(mp_up_diff8(mec_cells,:)),nanstd(mp_up_diff8(mec_cells,:))/sqrt(length(mec_cells)),'b');
hold on
h = errorbar(cwt_freqs,nanmean(mp_down_diff8(mec_cells,:)),nanstd(mp_down_diff8(mec_cells,:))/sqrt(length(mec_cells)),'r');
% h = errorbar(cwt_freqs,nanmean(mp_up_diff8(lec_cells,:)),nanstd(mp_up_diff8(lec_cells,:))/sqrt(length(lec_cells)),'k');
% h = errorbar(cwt_freqs,nanmean(mp_down_diff8(lec_cells,:)),nanstd(mp_down_diff8(lec_cells,:))/sqrt(length(lec_cells)),'g');
set(gca,'xscale','log')
xlim([2 100])
set(gca,'FontName','Arial','Fontsize',12)
xlabel('Frequency (Hz)','Fontsize',12,'fontname','Arial')
ylabel('Amplitude Difference (z)','Fontsize',12,'fontname','Arial')
legend('UP state','DOWN state')

%%
% figure
% h = errorbar(cwt_freqs,nanmean(lf8_up8_spec(mec_cells,:)),nanstd(lf8_up8_spec(mec_cells,:))/sqrt(length(mec_cells)),'b');
% hold on
% h = errorbar(cwt_freqs,nanmean(lf8_down8_spec(mec_cells,:)),nanstd(lf8_down8_spec(mec_cells,:))/sqrt(length(mec_cells)),'r');
% h = errorbar(cwt_freqs,nanmean(lf8_up_up8_spec(mec_cells,:)),nanstd(lf8_up_up8_spec(mec_cells,:))/sqrt(length(mec_cells)),'k');
% h = errorbar(cwt_freqs,nanmean(lf8_up_down8_spec(mec_cells,:)),nanstd(lf8_up_down8_spec(mec_cells,:))/sqrt(length(mec_cells)),'g');
% h = errorbar(cwt_freqs,nanmean(lf8_down_up8_spec(mec_cells,:)),nanstd(lf8_down_up8_spec(mec_cells,:))/sqrt(length(mec_cells)),'k');
% h = errorbar(cwt_freqs,nanmean(lf8_down_down8_spec(mec_cells,:)),nanstd(lf8_down_down8_spec(mec_cells,:))/sqrt(length(mec_cells)),'g');
% set(gca,'xscale','log')
% xlim([2 100])
% set(gca,'FontName','Arial','Fontsize',12)
% xlabel('Frequency (Hz)','Fontsize',12,'fontname','Arial')
% ylabel('Relative Amplitude (AU)','Fontsize',12,'fontname','Arial')
% 
% figure
% h = errorbar(cwt_freqs,nanmean(lf8_up8_spec(lec_cells,:)),nanstd(lf8_up8_spec(lec_cells,:))/sqrt(length(lec_cells)),'b');
% hold on
% h = errorbar(cwt_freqs,nanmean(lf8_down8_spec(lec_cells,:)),nanstd(lf8_down8_spec(lec_cells,:))/sqrt(length(lec_cells)),'r');
% h = errorbar(cwt_freqs,nanmean(lf8_up_up8_spec(lec_cells,:)),nanstd(lf8_up_up8_spec(lec_cells,:))/sqrt(length(lec_cells)),'k');
% h = errorbar(cwt_freqs,nanmean(lf8_up_down8_spec(lec_cells,:)),nanstd(lf8_up_down8_spec(lec_cells,:))/sqrt(length(lec_cells)),'g');
% h = errorbar(cwt_freqs,nanmean(lf8_down_up8_spec(lec_cells,:)),nanstd(lf8_down_up8_spec(lec_cells,:))/sqrt(length(lec_cells)),'k');
% h = errorbar(cwt_freqs,nanmean(lf8_down_down8_spec(lec_cells,:)),nanstd(lf8_down_down8_spec(mec_cells,:))/sqrt(length(lec_cells)),'g');
% set(gca,'xscale','log')
% xlim([2 100])
% set(gca,'FontName','Arial','Fontsize',12)
% xlabel('Frequency (Hz)','Fontsize',12,'fontname','Arial')
% ylabel('Relative Amplitude (AU)','Fontsize',12,'fontname','Arial')

lf8_up8_diff = lf8_up_up8_spec - lf8_down_up8_spec;
lf8_down8_diff = lf8_up_down8_spec - lf8_down_down8_spec;
lf8_up8_diff = lf8_up_up8_spec - lf8_down_up8_spec;
lf8_down8_diff = lf8_up_down8_spec - lf8_down_down8_spec;

figure
h = errorbar(cwt_freqs,nanmean(lf8_up8_diff(mec_cells,:)),nanstd(lf8_up8_diff(mec_cells,:))/sqrt(length(mec_cells)),'b');
hold on
h = errorbar(cwt_freqs,nanmean(lf8_down8_diff(mec_cells,:)),nanstd(lf8_down8_diff(mec_cells,:))/sqrt(length(mec_cells)),'r');
% h = errorbar(cwt_freqs,nanmean(lf8_up8_diff(lec_cells,:)),nanstd(lf8_up8_diff(lec_cells,:))/sqrt(length(lec_cells)),'k');
% h = errorbar(cwt_freqs,nanmean(lf8_down8_diff(lec_cells,:)),nanstd(lf8_down8_diff(lec_cells,:))/sqrt(length(lec_cells)),'g');
set(gca,'xscale','log')
xlim([2 100])
set(gca,'FontName','Arial','Fontsize',12)
xlabel('Frequency (Hz)','Fontsize',12,'fontname','Arial')
ylabel('Amplitude difference (z)','Fontsize',12,'fontname','Arial')
legend('UP state','DOWN state')

lf8_down_skipped_diff = lf8_down_skipped_spec - lf8_down_nonskipped_spec;
figure
h = errorbar(cwt_freqs,nanmean(lf8_down_skipped_diff(mec_cells,:)),nanstd(lf8_down_skipped_diff(mec_cells,:))/sqrt(length(mec_cells)),'b');
set(gca,'xscale','log')
xlim([1 100])
set(gca,'FontName','Arial','Fontsize',12)
xlabel('Frequency (Hz)','Fontsize',12,'fontname','Arial')
ylabel('Power difference (z)','Fontsize',12,'fontname','Arial')

lf8_down_long_diff = lf8_down_nonlong_spec-lf8_down_long_spec;
figure
h = errorbar(cwt_freqs,nanmean(lf8_down_long_diff(mec_cells,:)),nanstd(lf8_down_long_diff(mec_cells,:))/sqrt(length(mec_cells)),'b');
set(gca,'xscale','log')
xlim([1 100])
set(gca,'FontName','Arial','Fontsize',12)
xlabel('Frequency (Hz)','Fontsize',12,'fontname','Arial')
ylabel('Power difference (z)','Fontsize',12,'fontname','Arial')

figure
h = errorbar(nlx_amp_range,nanmean(lf8_down_skipped_amp(mec_cells,:)),nanstd(lf8_down_skipped_amp(mec_cells,:))./sqrt(sum(~isnan(lf8_down_skipped_amp))),'b');
errorbar_tick(h,.001,'units')
hold on
h = errorbar(nlx_amp_range,nanmean(lf8_down_nonskipped_amp(mec_cells,:)),nanstd(lf8_down_nonskipped_amp(mec_cells,:))./sqrt(sum(~isnan(lf8_down_nonskipped_amp))),'r');
errorbar_tick(h,.001,'units')
xlim([-2.5 1])
set(gca,'FontName','Arial','Fontsize',12)
xlabel('Amplitude (z)','Fontsize',12,'fontname','Arial')
ylabel('Probability','Fontsize',12,'fontname','Arial')
legend('Skipped','Non-Skipped')

%%
% figure
% h = errorbar(cwt_freqs,nanmean(lf3_up8_spec(mec_cells,:)),nanstd(lf3_up8_spec(mec_cells,:))/sqrt(length(mec_cells)),'b');
% hold on
% h = errorbar(cwt_freqs,nanmean(lf3_down8_spec(mec_cells,:)),nanstd(lf3_up8_spec(mec_cells,:))/sqrt(length(mec_cells)),'r');
% h = errorbar(cwt_freqs,nanmean(lf3_up_up8_spec(mec_cells,:)),nanstd(lf3_up_up8_spec(mec_cells,:))/sqrt(length(mec_cells)),'k');
% h = errorbar(cwt_freqs,nanmean(lf3_up_down8_spec(mec_cells,:)),nanstd(lf3_up_down8_spec(mec_cells,:))/sqrt(length(mec_cells)),'g');
% h = errorbar(cwt_freqs,nanmean(lf3_down_up8_spec(mec_cells,:)),nanstd(lf3_down_up8_spec(mec_cells,:))/sqrt(length(mec_cells)),'k');
% h = errorbar(cwt_freqs,nanmean(lf3_down_down8_spec(mec_cells,:)),nanstd(lf3_down_down8_spec(mec_cells,:))/sqrt(length(mec_cells)),'g');
% set(gca,'xscale','log')
% xlim([2 100])
% set(gca,'FontName','Arial','Fontsize',12)
% xlabel('Frequency (Hz)','Fontsize',12,'fontname','Arial')
% ylabel('Relative Amplitude (AU)','Fontsize',12,'fontname','Arial')
% 
% figure
% h = errorbar(cwt_freqs,nanmean(lf3_up8_spec(lec_cells,:)),nanstd(lf3_up8_spec(lec_cells,:))/sqrt(length(lec_cells)),'b');
% hold on
% h = errorbar(cwt_freqs,nanmean(lf3_down8_spec(lec_cells,:)),nanstd(lf3_up8_spec(lec_cells,:))/sqrt(length(lec_cells)),'r');
% h = errorbar(cwt_freqs,nanmean(lf3_up_up8_spec(lec_cells,:)),nanstd(lf3_up_up8_spec(lec_cells,:))/sqrt(length(lec_cells)),'k');
% h = errorbar(cwt_freqs,nanmean(lf3_up_down8_spec(lec_cells,:)),nanstd(lf3_up_down8_spec(lec_cells,:))/sqrt(length(lec_cells)),'g');
% h = errorbar(cwt_freqs,nanmean(lf3_down_up8_spec(lec_cells,:)),nanstd(lf3_down_up8_spec(lec_cells,:))/sqrt(length(lec_cells)),'k');
% h = errorbar(cwt_freqs,nanmean(lf3_down_down8_spec(lec_cells,:)),nanstd(lf3_down_down8_spec(mec_cells,:))/sqrt(length(lec_cells)),'g');
% set(gca,'xscale','log')
% xlim([2 100])
% set(gca,'FontName','Arial','Fontsize',12)
% xlabel('Frequency (Hz)','Fontsize',12,'fontname','Arial')
% ylabel('Relative Amplitude (AU)','Fontsize',12,'fontname','Arial')

lf3_up8_diff = lf3_up_up8_spec - lf3_down_up8_spec;
lf3_down8_diff = lf3_up_down8_spec - lf3_down_down8_spec;
lf3_up_diff8 = lf3_up_up8_spec - lf3_up_down8_spec;
lf3_down_diff8 = lf3_down_up8_spec - lf3_down_down8_spec;

figure
h = errorbar(cwt_freqs,nanmean(lf3_up8_diff(mec_cells,:)),nanstd(lf3_up8_diff(mec_cells,:))/sqrt(length(mec_cells)),'b');
hold on
h = errorbar(cwt_freqs,nanmean(lf3_down8_diff(mec_cells,:)),nanstd(lf3_down8_diff(mec_cells,:))/sqrt(length(mec_cells)),'r');
h = errorbar(cwt_freqs,nanmean(lf3_up_diff8(mec_cells,:)),nanstd(lf3_up_diff8(mec_cells,:))/sqrt(length(mec_cells)),'k');
h = errorbar(cwt_freqs,nanmean(lf3_down_diff8(mec_cells,:)),nanstd(lf3_down_diff8(mec_cells,:))/sqrt(length(mec_cells)),'g');
set(gca,'xscale','log')
xlim([2 100])
set(gca,'FontName','Arial','Fontsize',12)
xlabel('Frequency (Hz)','Fontsize',12,'fontname','Arial')
ylabel('Amplitude Difference (z)','Fontsize',12,'fontname','Arial')
legend('Cortial UP state','Cortical DOWN State','MP UP State','MP DOWN state')