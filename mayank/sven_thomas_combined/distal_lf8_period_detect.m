clear all
cd C:\WC_Germany\sven_thomas_combined
% load C:\WC_Germany\sven_thomas_combined\combined_dir.mat

load .\distal_dir
distal_dir = distal_dir(distal_usable);
mctx_lfp = ctx_lfp(distal_usable);
mhpc_lfp = hpc_lfp(distal_usable);
mhpc_mua = hpc_mua(distal_usable);

load ./distal_lec_dir.mat
distal_dir = [distal_dir distal_lec(usable_distal_lec)];
ctx_lfp = [mctx_lfp ctx_lfp(usable_distal_lec)];
hpc_lfp = [mhpc_lfp hpc_lfp(usable_distal_lec)];
hpc_mua = [mhpc_mua hpc_mua(usable_distal_lec)];
distal_mec = (1:length(distal_usable));
distal_lec = (length(distal_usable)+1):(length(distal_usable)+length(usable_distal_lec));
%%
dsf = 8;
Fsd = 2016/dsf;
niqf = 2016/2;
lcf = 0.05/niqf;
hcf = 2/niqf;
[b,a] = butter(2,[lcf hcf]);

%     lf8_up_fract = 0.4;
%     lf8_down_phase = 144;
up_fract = 0.5;
down_phase = 180;

for d = 1:length(distal_dir)
    
    fprintf('session %d\n',d)
    disp(distal_dir{d})
    cdir = distal_dir{d};
    cd(cdir)
    
    load ./used_data lf8 lf7 lf6 lf5
    if ctx_lfp(d) == 7
        lf8 = lf7;
    elseif ctx_lfp(d) == 6
        lf8 = lf6;
    elseif ctx_lfp(d) == 5
        lf8 = lf5;
    elseif ctx_lfp(d) ~= 8
        error('bad lfp num');
    end
    
    %     if ismember(d,new_data_inds)
    %         load ./pa_hsmm_state_seq8_newmec
    %         load ./pa_hsmm_state_seq_newmec
    %     else
    %         load ./pa_hsmm_state_seq_new2
    %         load ./pa_hsmm_state_seq8_new2
    %     end
    
    % load ./pa_hsmm_state_seq_combined_fin.mat
    % load ./pa_hsmm_state_seq7_combined_fin.mat
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


