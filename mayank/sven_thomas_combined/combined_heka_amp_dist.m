clear all

load C:\WC_Germany\sven_thomas_combined\combined_dir_nd
addpath('C:\WC_Germany\parietal_cortical_2010\')
addpath('C:\WC_Germany\hsmm_state_detection\')
addpath('C:\WC_Germany\persistent_2010\')
addpath('C:\Code\smoothing\software\')
addpath('C:\Code\general\')
l3mec_np(ismember(l3mec_np,[64])) = [];
uset = sort([l3mec l3lec l3mec_np l3lec_np]);
all_cells = 1:length(combined_dir);
l3mec = find(ismember(all_cells(uset),l3mec));
l3lec = find(ismember(all_cells(uset),l3lec));
l3mec_np = find(ismember(all_cells(uset),l3mec_np));
l3lec_np = find(ismember(all_cells(uset),l3lec_np));
combined_dir = combined_dir(uset);
% hpc_mua = hpc_mua(uset);
% hpc_lfp = hpc_lfp(uset);
ctx_lfp = ctx_lfp(uset);

ac_Fs = 2016;
ac_dsf = 8;
ac_Fsd = 2016/ac_dsf;
ac_niqf = 2016/2;
[b,a] = butter(2,[0.05/ac_niqf 10/ac_niqf]);
[b2,a2] = butter(2,[0.05/ac_niqf 40/ac_niqf]);

hf_lcf = 15;
hf_hcf = 100;
[b3,a3] = butter(2,[hf_lcf hf_hcf]/ac_niqf);

dc_target_fs = 200;

heka_amp_range = linspace(-90,-20,500);
nlx_amp_range = linspace(-4,4,500);
%%
for d = 1:length(combined_dir)
% for d = 30
    cd(combined_dir{d});
    d
    load ./used_data lf8 lf7 wcv_minus_spike
    if ctx_lfp(d) == 7
        lf8 = lf7;
    end
    lf8_d = zscore(downsample(filtfilt(b2,a2,lf8),ac_dsf));
    
    load ./aligned_heka
    load ./heka_aligned_spksub
    
    wcv_d = zscore(downsample(filtfilt(b2,a2,wcv_minus_spike),ac_dsf));
    
    dc_Fs = 1/(median(diff(dc_time)));
    dc_dsf = round(dc_Fs/dc_target_fs);
    dc_data = downsample(dc_data,dc_dsf);
    dc_data_spksub = downsample(dc_data_spksub,dc_dsf);
    
    if mean(abs(dc_data)) < 1
        dc_data = dc_data*100;
    end
    
    dc_time = downsample(dc_time,dc_dsf);
    ac_time = downsample(ac_time,ac_dsf);
    
    mp_state_var = nan(size(wcv_d));
    mp_state_var_sig = nan(size(wcv_d));
    
    %     load ./pa_hsmm_state_seq_combined_fin
    load ./pa_hsmm_state_seq_combined_fin_nd
    mp_state_seq = hsmm_bbstate_seq;
    
    dc_data = dc_data(:);
    dc_data_spksub = dc_data_spksub(:);
    dc_upvals = [];
    dc_downvals = [];
    dc_upvals_spksub = [];
    dc_downvals_spksub = [];
    dc_upvals_spksub_sig = [];
    dc_downvals_spksub_sig = [];
    
    [new_seg_inds] = resample_uds_seg_inds_v2(hmm.UDS_segs,hmm.Fs,ac_Fsd,mp_state_seq);
    [up_trans_inds,down_trans_inds] = compute_state_transitions_seg(new_seg_inds,mp_state_seq);
    
    mp_upstate_heka_meanamp{d} = nan(size(up_trans_inds));
    mp_downstate_heka_meanamp{d} = nan(size(up_trans_inds));
    cnt = 1;
    for i = 1:length(mp_state_seq)
        up_trans = find(mp_state_seq{i}(1:end-1)==1 & mp_state_seq{i}(2:end)==2);
        down_trans = find(mp_state_seq{i}(1:end-1)==2 & mp_state_seq{i}(2:end)==1);
        up_trans(up_trans > down_trans(end)) = [];
        down_trans(down_trans < up_trans(1)) = [];
        up_trans = up_trans + new_seg_inds(i,1)-1;
        down_trans = down_trans + new_seg_inds(i,1)-1;
        up_trans = round(up_trans/(ac_dsf/8));
        down_trans = round(down_trans/(ac_dsf/8));
        for j = 1:length(up_trans)
            cur_vals = find(dc_time > ac_time(up_trans(j)) & dc_time < ac_time(down_trans(j)));
            dc_upvals = [dc_upvals; dc_data(cur_vals)];
            dc_upvals_spksub = [dc_upvals_spksub; dc_data_spksub(cur_vals)];
            mp_upstate_heka_meanamp{d}(cnt) = nanmean(dc_data_spksub(cur_vals));
            if j < length(up_trans)
                cur_vals = find(dc_time > ac_time(down_trans(j)) & dc_time < ac_time(up_trans(j+1)));
                dc_downvals = [dc_downvals; dc_data(cur_vals)];
                dc_downvals_spksub = [dc_downvals_spksub; dc_data_spksub(cur_vals)];
                mp_downstate_heka_meanamp{d}(cnt) = nanmean(dc_data_spksub(cur_vals));
            end
            
            cnt = cnt + 1;
        end
        mp_state_var(new_seg_inds(i,1):new_seg_inds(i,2)) = mp_state_seq{i};
    end
    
    %manually get rid of chunk where state detection clearly failing
    if d == 5
        dc_downvals_spksub(3000:13000) = [];
    end
    
    %     beg_states = find(up_trans_inds < up_trans_inds(1) + 50*ac_Fsd);
    %     upstate_mean_spksub_beg(d) = nanmean(mp_upstate_heka_meanamp{d}(beg_states));
    %     downstate_mean_spksub_beg(d) = nanmean(mp_downstate_heka_meanamp{d}(beg_states));
    
    use_up_pts = find(~isnan(dc_upvals_spksub));
    use_down_pts = find(~isnan(dc_downvals_spksub));
    beg_down_data = 1:round(length(dc_downvals_spksub)/20);
    beg_up_data = 1:round(length(dc_upvals_spksub)/20);
    upstate_mean_spksub_beg(d) = nanmean(dc_upvals_spksub(beg_up_data));
    downstate_mean_spksub_beg(d) = nanmean(dc_downvals_spksub(beg_down_data));
    upstate_mode_spksub_beg(d) = mode(dc_upvals_spksub(beg_up_data));
    downstate_mode_spksub_beg(d) = mode(dc_downvals_spksub(beg_down_data));
    upstate_mode_spksub_beg(d) = median(dc_upvals_spksub(beg_up_data));
    downstate_mode_spksub_beg(d) = median(dc_downvals_spksub(beg_down_data));
    upstate_mean(d) = nanmean(dc_upvals);
    downstate_mean(d) = nanmean(dc_downvals);
    upstate_var(d) = nanvar(dc_upvals);
    downstate_var(d) = nanvar(dc_downvals);
    upstate_mode_spksub(d) = mode(dc_upvals_spksub(use_up_pts));
    downstate_mode_spksub(d) = mode(dc_downvals_spksub(use_down_pts));
    upstate_mean_spksub(d) = nanmean(dc_upvals_spksub);
    downstate_mean_spksub(d) = nanmean(dc_downvals_spksub);
    upstate_median_spksub(d) = nanmean(dc_upvals_spksub);
    downstate_median_spksub(d) = nanmean(dc_downvals_spksub);
    upstate_var_spksub(d) = nanvar(dc_upvals_spksub);
    downstate_var_spksub(d) = nanvar(dc_downvals_spksub);
    
    upstate_heka_dist(d,:) = ksdensity(dc_upvals_spksub,heka_amp_range);
    downstate_heka_dist(d,:) = ksdensity(dc_downvals_spksub,heka_amp_range);
    ov_heka_dist(d,:) = ksdensity([dc_upvals(:); dc_downvals(:)],heka_amp_range);
    
    upstate_nlx_dist(d,:) = ksdensity(wcv_d(mp_state_var==2),nlx_amp_range);
    downstate_nlx_dist(d,:) = ksdensity(wcv_d(mp_state_var==1),nlx_amp_range);
    ov_nlx_dist(d,:) = ksdensity(wcv_d(~isnan(mp_state_var)),nlx_amp_range);
    
    clear hsmm* wcv* lf*
end

cd C:\WC_Germany\sven_thomas_combined\
save combined_heka_UDS_data_allcells_fin_nd_np upstate* downstate* ov_* *amp_range mp_*state*

%%
% Y = [

