clear all

load F:\WC_Germany\overall_EC\overall_EC_dir.mat
addpath('F:\WC_Germany\parietal_cortical_2010\')
addpath('F:\WC_Germany\hsmm_state_detection\')
addpath('F:\WC_Germany\persistent_2010\')
addpath('F:\Code\smoothing\software\')
addpath('F:\Code\general\')
addpath('F:\WC_Germany\persistent_9_27_2010\')

used_data = [l3mec l3lec];
sess_data = sess_data(used_data);

Fs = 2016;
dsf = 8;
Fsd = 2016/dsf;
niqf = Fs/2;
[b,a] = butter(2,[0.05/niqf 20/niqf]);

min_rsquared = 0.6;

for d = 7:10
    cdir = sess_data(d).directory;
    cdir(1) = 'F';
    disp(sprintf('session %d',d))
    cd(cdir);
    s_name = strcat(sess_data(d).region,'_l',sess_data(d).layer,'_',sess_data(d).name);
    
%     if ~exist('pa_sig_fit_data.mat','file')
        
        load ./used_data lf8 wcv_minus_spike
        wcv_d = zscore(downsample(filtfilt(b,a,wcv_minus_spike),dsf));
        lf8_d = zscore(downsample(filtfilt(b,a,lf8),dsf));
        t_axis = (1:length(wcv_d))/Fsd;
        
        load ./pa_hsmm_state_seq
        mp_state_seq = hsmm_bbstate_seq;
        
        [new_seg_inds] = resample_uds_seg_inds(hmm.UDS_segs,hmm.Fs,252,length(wcv_d));
        [up_trans_inds,down_trans_inds] = compute_state_transitions_seg(new_seg_inds,mp_state_seq);
        
        [ruamp,rushift,rutau,rusquared,tu_50,tu_90,tu_10,ufit_data] = ...
            get_sigmoid_fit_ut_pa(up_trans_inds,down_trans_inds,wcv_d,t_axis,Fsd);
        used_ufits = find(~isnan(rusquared) & rusquared > min_rsquared);
        
        [rdamp,rdshift,rdtau,rdsquared,td_50,td_90,td_10,dfit_data] = ...
            get_sigmoid_fit_dt_pa(up_trans_inds,down_trans_inds,wcv_d,t_axis,Fsd);
        used_dfits = find(~isnan(rdsquared) & rdsquared > min_rsquared);
        
        save pa_sig_fit_data ru* tu_* rd* td_*
        clear ru* tu_* rd* td_*
        
        load ./pa_hsmm_state_seq8
        lfp_state_seq = hsmm_bbstate_seq8;
        
        [new_seg_inds] = resample_uds_seg_inds(hmm8.UDS_segs,hmm8.Fs,252,length(lf8_d));
        [up_trans_inds8,down_trans_inds8] = compute_state_transitions_seg(new_seg_inds,lfp_state_seq);
        
        [ruamp8,rushift8,rutau8,rusquared8,tu8_50,tu8_90,tu8_10,ufit_data8] = ...
            get_sigmoid_fit_ut_pa(up_trans_inds8,down_trans_inds8,lf8_d,t_axis,Fsd);
        used_ufits8 = find(~isnan(rusquared8) & rusquared8 > min_rsquared);
        
        [rdamp8,rdshift8,rdtau8,rdsquared8,td8_50,td8_90,td8_10,dfit_data8] = ...
            get_sigmoid_fit_dt_pa(up_trans_inds8,down_trans_inds8,lf8_d,t_axis,Fsd);
        used_dfits8 = find(~isnan(rdsquared8) & rdsquared8 > min_rsquared);
        
        save pa_sig_fit_data8 ru* tu8_* rd* td8_*
        clear ru* tu8_* rd* td8_*
        
        %     figure
        %     plot(t_axis,wcv_d), hold on
        % %     plot(t_axis(t_90(used_fits)),wcv_d(t_90(used_fits)),'ro')
        % %     plot(t_axis(t_10(used_fits)),wcv_d(t_10(used_fits)),'go')
        %     plot(t_axis,ufit_data,'r','linewidth',2)
        %     plot(t_axis,dfit_data,'k','linewidth',2)
        %     pause
        %     close all
%     end
end

cd F:\WC_Germany\persistent_9_27_2010\\
