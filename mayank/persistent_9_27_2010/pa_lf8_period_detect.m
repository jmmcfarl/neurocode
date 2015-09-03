clear all

load G:\WC_Germany\overall_EC\overall_EC_dir
drive_letter = 'G';
used_data = [l3mec_p l3lec_p];
sess_data = sess_data(used_data);

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

for d = 1:length(sess_data)
    
    fprintf('session %d\n',d)
    disp(sess_data(d).name)
    cdir = sess_data(d).directory;
    cdir(1) = drive_letter;
    cd(cdir)
    
%     if ~exist('./pa_lf8_period_data.mat','file')
        load ./used_data lf8
%         load ./pa_hsmm_state_seq8_sm500
% %         lfp_state_seq = hmm_bbstate_seq8;
        load ./pa_hsmm_state_seq_new2
        lfp_state_seq = hmm_bbstate_seq;
        hmm8 = hmm;
        
        lf8_f = filtfilt(b,a,lf8);
        lf8_d = zscore(downsample(lf8_f,dsf));
             
        [new_seg_inds] = resample_uds_seg_inds(hmm8.UDS_segs,hmm8.Fs,Fsd,length(lf8_d));
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
        save pa_lf8_period_data_new2_reverse lf8_period* lf8_uds_amp
        clear lf8*
%     end
end


