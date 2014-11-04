function [is_bad] = all_EC_compute_ctx_period()

%%
dsf = 8;
Fsd = 2016/dsf;
niqf = 2016/2;
lcf = 0.05/niqf;
hcf = 2/niqf;
[b,a] = butter(2,[lcf hcf]);

up_fract = 0.5;
down_phase = 180;

load ./used_data lf7

% load ./all_combined_mp_uds.mat
% load ./all_combined_lf7_uds.mat
load ./pa_hsmm_state_seq_combined_fin_nd.mat
load ./pa_hsmm_state_seq7_combined_fin_nd.mat

hmm8 = hmm7;
lfp_state_seq = hmm_bbstate_seq7;

lf8_f = filtfilt(b,a,lf7);
lf8_d = zscore(downsample(lf8_f,dsf));

if isempty(hmm8)
    is_bad = 1;
else
    is_bad = 0;
    [new_seg_inds] = resample_uds_seg_inds_v2(hmm8.UDS_segs,hmm8.Fs,Fsd,lfp_state_seq);
    seg_durs = new_seg_inds(:,2)-new_seg_inds(:,1)+1;
    
    
    for ns = 1:hmm8.Nsegs
        cur_up_trans = find(lfp_state_seq{ns}(1:end-1) == 1 & lfp_state_seq{ns}(2:end) == 2);
        cur_down_trans = find(lfp_state_seq{ns}(1:end-1) == 2 & lfp_state_seq{ns}(2:end) == 1);
        cur_up_trans(cur_up_trans > cur_down_trans(end)) = [];
        cur_down_trans(cur_down_trans < cur_up_trans(1)) = [];
        
        lf8_period_f{ns} = nan(seg_durs(ns),1);
        lf8_period_p{ns} = nan(seg_durs(ns),1);
        lf8_uds_amp{ns} = nan(length(cur_up_trans)-1,1);
        for i = 1:length(cur_up_trans)-1
            cur_up_times = cur_up_trans(i):cur_down_trans(i);
            cur_down_times = cur_down_trans(i):cur_up_trans(i+1);
            max_up = max(lf8_d(cur_up_times));
            min_down = min(lf8_d(cur_down_times));
            lf8_uds_amp{ns}(i) = max_up-min_down;
            
            period_samps_up = cur_down_trans(i)-cur_up_trans(i);
            period_samps_down = cur_up_trans(i+1)-cur_down_trans(i);
            
            lf8_period_f{ns}(cur_up_trans(i)+1:cur_down_trans(i)) = ...
                i-1+linspace(1,period_samps_up,period_samps_up)/period_samps_up*up_fract;
            lf8_period_f{ns}(cur_down_trans(i)+1:cur_up_trans(i+1)) = ...
                i-1+up_fract+linspace(1,period_samps_down,period_samps_down)/period_samps_down*(1-up_fract);
            
            lf8_period_p{ns}(cur_up_trans(i)+1:cur_down_trans(i)) = ...
                linspace(1,period_samps_up,period_samps_up)/period_samps_up*down_phase;
            lf8_period_p{ns}(cur_down_trans(i)+1:cur_up_trans(i+1)) = ...
                down_phase+linspace(1,period_samps_down,period_samps_down)/period_samps_down*(360-down_phase);
            
        end
    end
%     save allEC_ctx_period_data lf8_period* lf8_uds_amp
    save allEC_ctx_period_data_hsmm lf8_period* lf8_uds_amp
    clear lf8* hmm8
end


