clear all

load G:\WC_Germany\overall_EC\overall_EC_dir.mat
addpath('G:\WC_Germany\parietal_cortical_2010\')
addpath('G:\WC_Germany\hsmm_state_detection\')
addpath('G:\WC_Germany\persistent_2010\')
addpath('G:\Code\smoothing\software\')
addpath('G:\Code\general\')

used_data = [l3mec_p l3lec_p];
sess_data = sess_data(used_data);

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

min_sig_rsquared = 0.6;

for d = 1:length(sess_data)
    cdir = sess_data(d).directory;
    cdir(1) = 'G';
    disp(sprintf('session %d',d))
    cd(cdir);
    s_name = strcat(sess_data(d).region,'_l',sess_data(d).layer,'_',sess_data(d).name);
    
    load ./used_data lf8 wcv_minus_spike
    lf8_d = zscore(downsample(filtfilt(b2,a2,lf8),ac_dsf));
  
    if ~isempty(sess_data(d).heka_dir)
        
        load ./aligned_heka
        load ./heka_aligned_spksub
        load ./pa_sig_fit_data
        load ./pa_sig_fit_data8
        
        wcv_d = zscore(downsample(filtfilt(b2,a2,wcv_minus_spike),ac_dsf));
        wcv_hf = downsample(filtfilt(b3,a3,wcv_minus_spike),ac_dsf);
        
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
        
        load ./pa_hsmm_state_seq_sm500
        mp_state_seq = hsmm_bbstate_seq;
        dc_data = dc_data(:);
        dc_data_spksub = dc_data_spksub(:);
        dc_upvals = [];
        dc_downvals = [];
        dc_upvals_spksub = [];
        dc_downvals_spksub = [];
        dc_upvals_spksub_sig = [];
        dc_downvals_spksub_sig = [];
        
        [new_seg_inds] = resample_uds_seg_inds(hmm.UDS_segs,hmm.Fs,ac_Fsd,length(wcv_d));
        [up_trans_inds,down_trans_inds] = compute_state_transitions_seg(new_seg_inds,mp_state_seq);
        
        mp_upstate_heka_meanamp{d} = nan(size(up_trans_inds));
        mp_downstate_heka_meanamp{d} = nan(size(up_trans_inds));
        mp_upstate_nlx_stdamp{d} = nan(size(up_trans_inds));
        mp_downstate_nlx_stdamp{d} = nan(size(up_trans_inds));
%         mp_upstate_heka_meanamp_sig{d} = nan(size(up_trans_inds));
%         mp_downstate_heka_meanamp_sig{d} = nan(size(up_trans_inds));
%         mp_upstate_nlx_stdamp_sig{d} = nan(size(up_trans_inds));
%         mp_downstate_nlx_stdamp_sig{d} = nan(size(up_trans_inds));
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
                mp_upstate_nlx_stdamp{d}(cnt) = nanstd(wcv_hf(up_trans(j):down_trans(j)));
%                 is_sig_good = rusquared(j) > min_sig_rsquared & rdsquared(j) > min_sig_rsquared;
%                 if is_sig_good
%                     cur_sig_vals = find(dc_time > ac_time(tu_90(j)) & dc_time < ac_time(td_10(j)));
%                     mp_upstate_heka_meanamp_sig{d}(cnt) = nanmean(dc_data_spksub(cur_sig_vals));
%                     mp_upstate_nlx_stdamp_sig{d}(cnt) = nanstd(wcv_hf(tu_90(j):td_10(j)));
%                     dc_upvals_spksub_sig = [dc_upvals_spksub_sig; dc_data_spksub(cur_sig_vals)];
%                     mp_state_var_sig(tu_90(j):td_10(j)) = 2;
%                 end
                if j < length(up_trans)
                    cur_vals = find(dc_time > ac_time(down_trans(j)) & dc_time < ac_time(up_trans(j+1)));
                    dc_downvals = [dc_downvals; dc_data(cur_vals)];
                    dc_downvals_spksub = [dc_downvals_spksub; dc_data_spksub(cur_vals)];
                    mp_downstate_heka_meanamp{d}(cnt) = nanmean(dc_data_spksub(cur_vals));
                    mp_downstate_heka_stdamp{d}(cnt) = nanstd(wcv_hf(down_trans(j):up_trans(j+1)));
%                     is_sig_good = rdsquared(j) > min_sig_rsquared & rusquared(j+1) > min_sig_rsquared;
%                     if is_sig_good
%                         cur_sig_vals = find(dc_time > ac_time(td_90(j)) & dc_time < ac_time(tu_10(j+1)));
%                         mp_downstate_heka_meanamp_sig{d}(cnt) = nanmean(dc_data_spksub(cur_sig_vals));
%                         mp_downstate_nlx_stdamp_sig{d}(cnt) = nanstd(wcv_hf(td_90(j):tu_10(j+1)));
%                         dc_downvals_spksub_sig = [dc_downvals_spksub_sig; dc_data_spksub(cur_sig_vals)];
%                         mp_state_var_sig(td_90(j):tu_10(j+1)) = 1;
%                     end
                end
                
                cnt = cnt + 1;
            end
            mp_state_var(new_seg_inds(i,1):new_seg_inds(i,2)) = mp_state_seq{i};
        end
        
        upstate_mean(d) = mean(dc_upvals);
        downstate_mean(d) = mean(dc_downvals);
        upstate_var(d) = var(dc_upvals);
        downstate_var(d) = var(dc_downvals);
        upstate_mean_spksub(d) = mean(dc_upvals_spksub);
        downstate_mean_spksub(d) = mean(dc_downvals_spksub);
%         upstate_mean_spksub_sig(d) = mean(dc_upvals_spksub_sig);
%         downstate_mean_spksub_sig(d) = mean(dc_downvals_spksub_sig);
        upstate_var_spksub(d) = var(dc_upvals_spksub);
        downstate_var_spksub(d) = var(dc_downvals_spksub);
%         upstate_var_spksub_sig(d) = var(dc_upvals_spksub_sig);
%         downstate_var_spksub_sig(d) = var(dc_downvals_spksub_sig);
        
        upstate_heka_dist(d,:) = ksdensity(dc_upvals,heka_amp_range);
        downstate_heka_dist(d,:) = ksdensity(dc_downvals,heka_amp_range);
%         upstate_heka_dist_spksub_sig(d,:) = ksdensity(dc_upvals_spksub_sig,heka_amp_range);
%         downstate_heka_dist_spksub_sig(d,:) = ksdensity(dc_downvals_spksub_sig,heka_amp_range);
        ov_heka_dist(d,:) = ksdensity([dc_upvals(:); dc_downvals(:)],heka_amp_range);
        
        upstate_nlx_dist(d,:) = ksdensity(wcv_d(mp_state_var==2),nlx_amp_range);
        downstate_nlx_dist(d,:) = ksdensity(wcv_d(mp_state_var==1),nlx_amp_range);
        ov_nlx_dist(d,:) = ksdensity(wcv_d(~isnan(mp_state_var)),nlx_amp_range);
        
        upstate_var_hfnlx(d) = var(wcv_hf(mp_state_var==2));
        downstate_var_hfnlx(d) = var(wcv_hf(mp_state_var==1));
%         upstate_var_hfnlx_sig(d) = var(wcv_hf(mp_state_var_sig==2));
%         downstate_var_hfnlx_sig(d) = var(wcv_hf(mp_state_var_sig==1));
        
    end
    
    %%
    load ./pa_hsmm_state_seq8_sm500
    lfp_state_var = nan(size(lf8_d));
%     lfp_state_var_sig = nan(size(lf8_d));
    lfp_state_seq = hsmm_bbstate_seq8;
    [new_seg_inds] = resample_uds_seg_inds(hmm8.UDS_segs,hmm8.Fs,252,length(lf8_d));
    [up_trans_inds8,down_trans_inds8] = compute_state_transitions_seg(new_seg_inds,lfp_state_seq);
    cnt = 1;
    for i = 1:length(lfp_state_seq)
        lfp_state_var(new_seg_inds(i,1):new_seg_inds(i,2)) = lfp_state_seq{i};
        up_trans = find(lfp_state_seq{i}(1:end-1)==1 & lfp_state_seq{i}(2:end)==2);
        down_trans = find(lfp_state_seq{i}(1:end-1)==2 & lfp_state_seq{i}(2:end)==1);
        up_trans(up_trans > down_trans(end)) = [];
        down_trans(down_trans < up_trans(1)) = [];
        up_trans = up_trans + new_seg_inds(i,1)-1;
        down_trans = down_trans + new_seg_inds(i,1)-1;
        up_trans = round(up_trans/(ac_dsf/8));
        down_trans = round(down_trans/(ac_dsf/8));
%         for j = 1:length(up_trans)
%             is_sig_good = rusquared8(j) > min_sig_rsquared & rdsquared8(j) > min_sig_rsquared;
%             if is_sig_good
%                 lfp_state_var_sig(tu8_90(j):td8_10(j)) = 2;
%             end
%             if j < length(up_trans)
%                 is_sig_good = rdsquared8(j) > min_sig_rsquared & rusquared8(j+1) > min_sig_rsquared;
%                 if is_sig_good
%                     lfp_state_var_sig(td8_90(j):tu8_10(j+1)) = 1;
%                 end
%             end
%         end
    end
    
    upstate_nlx_dist8(d,:) = ksdensity(lf8_d(lfp_state_var==2),nlx_amp_range);
    downstate_nlx_dist8(d,:) = ksdensity(lf8_d(lfp_state_var==1),nlx_amp_range);
%     upstate_nlx_dist8_sig(d,:) = ksdensity(lf8_d(lfp_state_var_sig==2),nlx_amp_range);
%     downstate_nlx_dist8_sig(d,:) = ksdensity(lf8_d(lfp_state_var_sig==1),nlx_amp_range);
    ov_nlx_dist8(d,:) = ksdensity(lf8_d(~isnan(lfp_state_var)),nlx_amp_range);
    
    clear hsmm* wcv* lf*
    
end

cd G:\WC_Germany\persistent_9_27_2010\
save pa_heka_UDS_data_allcells_sm500 upstate* downstate* ov_* *amp_range mp_*state*