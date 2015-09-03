clear all

load G:\WC_Germany\overall_EC\overall_EC_dir.mat

used_data = [l3mec_p l3lec_p];
sess_data = sess_data(used_data);


dsf = 8;
Fsd = 2016/dsf;

for d = 1:length(sess_data)
    
    cdir = sess_data(d).directory;
    cdir(1) = 'G';
    disp(sprintf('session %d',d))
    cd(cdir);
    s_name = strcat(sess_data(d).region,'_l',sess_data(d).layer,'_',sess_data(d).name);
    
    load ./spike_time_jmm
    load ./used_data wcv
    
    wcv_d = downsample(wcv,dsf);
    spike_ids = round(spkid/dsf);
    ov_spike_log = zeros(size(wcv_d));
    ov_spike_log(spike_ids) = 1;

    load ./pa_hsmm_state_seq_sm500
    mp_state_seq = hsmm_bbstate_seq;
    load ./pa_hsmm_state_seq8_sm500
    lfp_state_seq = hsmm_bbstate_seq8;
    
    [new_seg_inds] = resample_uds_seg_inds(hmm.UDS_segs,hmm.Fs,Fsd,length(wcv_d));
    seg_durs = new_seg_inds(:,2)-new_seg_inds(:,1)+1;
    
    [up_trans_inds,down_trans_inds] = compute_state_transitions_seg(new_seg_inds,mp_state_seq);
    up_state_durs = (down_trans_inds-up_trans_inds)/Fsd;
    
    wcv_up_log = nan(size(wcv_d));
    lf8_up_log = nan(size(wcv_d));
    spike_log = nan(size(wcv_d));
    
    for ns = 1:hmm.Nsegs
        cur_seg = new_seg_inds(ns,1):new_seg_inds(ns,2);
        wcv_up_log(cur_seg) = logical(mp_state_seq{ns}-1);
        lf8_up_log(cur_seg) = logical(lfp_state_seq{ns}-1);
        spike_log(cur_seg) = ov_spike_log(cur_seg);
    end
    
    mp_upstate_rate{d} = nan(size(up_trans_inds));
    for j = 1:length(up_trans_inds)
        mp_upstate_rate{d}(j) = sum(spike_log(up_trans_inds(j):down_trans_inds(j)))/up_state_durs(j);
    end
    
    mp_up_rate(d) = nansum(spike_log(wcv_up_log==1))/nansum(wcv_up_log)*Fsd;
    mp_up_rate_lfup(d) = nansum(spike_log(wcv_up_log==1 & lf8_up_log==1))/nansum(wcv_up_log==1 & lf8_up_log==1)*Fsd;
    mp_up_rate_lfdown(d) = nansum(spike_log(wcv_up_log==1 & lf8_up_log==0))/nansum(wcv_up_log==1 & lf8_up_log==0)*Fsd;
    mp_rate(d) = nansum(spike_log)/length(spike_log)*Fsd;
    
    clear wcv wcv_d
       
end

mod_index = (mp_up_rate_lfup-mp_up_rate_lfdown)./(mp_up_rate_lfup+mp_up_rate_lfdown);

cd G:\WC_Germany\persistent_9_27_2010\
save spike_rate_data_sm500 mp_*

% 
mec_cells = 1:length(l3mec_p);
lec_cells = (mec_cells(end)+1):(mec_cells(end)+length(l3lec_p));

plot(mp_up_rate_lfup(mec_cells),mp_up_rate_lfdown(mec_cells),'.','markersize',14)
hold on
plot(mp_up_rate_lfup(lec_cells),mp_up_rate_lfdown(lec_cells),'g.','markersize',14)
line([0 12],[0 12],'Color','k')