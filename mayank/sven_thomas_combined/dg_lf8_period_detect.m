clear all
load C:\WC_Germany\sven_thomas_combined\dg_dir
%%
dsf = 8;
Fsd = 2016/dsf;
niqf = 2016/2;
lcf = 0.05/niqf;
hcf = 2/niqf;
[b,a] = butter(2,[lcf hcf]);

up_fract = 0.5;
down_phase = 180;

for d = 1:length(dg_dir)
    
    fprintf('session %d\n',d)
    disp(dg_dir{d})
    cdir = dg_dir{d};
    cd(cdir)
    
    load ./used_data lf8 lf7 
        lf8 = lf7;
    
    load ./pa_hsmm_state_seq_combined_fin_nd.mat
    load ./pa_hsmm_state_seq7_combined_fin_nd.mat
    hsmm8 = hsmm7;
    lfp_state_seq = hsmm_bbstate_seq7;
    
%     %for reversing the roles of the MP and LFP
%     lfp_state_seq = hsmm_bbstate_seq;
    
    lf8_f = filtfilt(b,a,lf8);
    lf8_d = zscore(downsample(lf8_f,dsf));
    
    [new_seg_inds] = resample_uds_seg_inds_v2(hsmm8.UDS_segs,hsmm8.Fs,Fsd,lfp_state_seq);
    seg_durs = new_seg_inds(:,2)-new_seg_inds(:,1)+1;
    
    for ns = 1:hsmm8.Nsegs
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
        save combined_lf7_period_data_fin_nd lf8_period* lf8_uds_amp
%     save combined_mp_period_data_fin_nd lf8_period* lf8_uds_amp
    clear lf8*
end


