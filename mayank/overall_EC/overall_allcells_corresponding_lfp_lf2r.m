clear all

load F:\WC_Germany\overall_EC\overall_allcells_dir.mat

drive_letter = 'F';
dsf = 8;
raw_Fs = 2016;
Fsd = raw_Fs/dsf;
lf8_up_fract = 0.5;

for d = 1:length(sess_data)
    
    fprintf('session %d\n',d)
    cdir = sess_data(d).directory;
    cdir(1) = drive_letter;
    cd(cdir)
    
    load used_data lf2 lf5
    load desynch_times_lf2r
%     load desynch_times_lf3
    lf2 = lf2/sess_data(d).gains(2);
    lf5 = lf5/sess_data(d).gains(5);
    lf2_r = downsample(lf2-lf5,dsf);
%     lf3 = downsample(lf3,dsf)/sess_data(d).gains(3);
    total_time = length(lf2_r)/Fsd;
    desynch_time_2r = sum(desynch_times_lf2r(:,2)-desynch_times_lf2r(:,1));
%     desynch_time_3 = sum(desynch_times_lf3(:,2)-desynch_times_lf3(:,1));
    uds_time_2r = total_time - desynch_time_2r;
%     uds_time_3 = total_time - desynch_time_3;
    
    if max(isnan(lf2_r)) > 0
        bad_lf2r(d) = 1;
    else
        bad_lf2r(d) = 0;
    end
%     if max(isnan(lf3)) > 0
%         bad_lf3(d) = 1;
%     else
%         bad_lf3(d) = 0;
%     end
    
    load ec_hmm_state_seq8
    lfp_state_seq = hmm_bbstate_seq8;
    [new_seg_inds8] = resample_uds_seg_inds(hmm8.UDS_segs,hmm8.Fs,Fsd,length(lf2_r));
   
    load lf8_period_data
    is_desynch8 = ones(size(lf2_r));
    for ns = 1:hmm8.Nsegs
        is_desynch8(new_seg_inds8(ns,1):new_seg_inds8(ns,2)) = 0;
    end
    ov_lf8_period_f = nan(size(lf2_r));
    ov_lf8_period_p = nan(size(lf2_r));
    ov_up_trans8 = [];
    ov_down_trans8 = [];
    for ns = 1:hmm8.Nsegs
        cur_up_trans8 = new_seg_inds8(ns,1)-1+find(lfp_state_seq{ns}(1:end-1) == 1 & lfp_state_seq{ns}(2:end) == 2);
        cur_down_trans8 = new_seg_inds8(ns,1)-1+find(lfp_state_seq{ns}(1:end-1) == 2 & lfp_state_seq{ns}(2:end) == 1);
        cur_up_trans8(cur_up_trans8 > cur_down_trans8(end)) = [];
        cur_down_trans8(cur_down_trans8 < cur_up_trans8(1)) = [];
        ov_up_trans8 = [ov_up_trans8; cur_up_trans8];
        ov_down_trans8 = [ov_down_trans8; cur_down_trans8];
        ov_lf8_period_f(new_seg_inds8(ns,1):new_seg_inds8(ns,2)) = lf8_period_f{ns};
        ov_lf8_period_p(new_seg_inds8(ns,1):new_seg_inds8(ns,2)) = lf8_period_p{ns};
    end
    n_lfp_ups = length(ov_up_trans8);
   
    
    if bad_lf2r(d) == 0 && uds_time_2r > 250
        
        datalen = length(lf2_r);
        load ec_hmm_state_seq2r
        lf2r_state_seq = hmm_bbstate_seq2r;
        
        [new_seg_inds] = resample_uds_seg_inds(hmm2r.UDS_segs,hmm2r.Fs,Fsd,length(lf2_r));
        seg_durs = new_seg_inds(:,2)-new_seg_inds(:,1)+1;
        
        
        lf2r_up_inds{d} = [];
        lf2r_down_inds{d} = [];
        lf2r_updur{d} = [];
        lf2r_up_lag{d} = [];
        lf2r_down_lag{d} = [];
        lf2r_updur_lfpc{d} = [];
        corresp_lfp_updur{d} = [];
        corresp_lfp_downdur{d} = [];
        corresp_lfp_cycledur{d} = [];
        corresp_lfp_dutycycle{d} = [];
        corresp_lfp_alf2r{d} = [];
        corresp_lfp_ind{d} = [];
        lf2r_skipped_lfpupind{d} = [];
        lup_mdown_lag{d} = [];
                
        for ns = 1:hmm2r.Nsegs
            
            cur_up_trans = new_seg_inds(ns,1)-1+find(lf2r_state_seq{ns}(1:end-1) == 1 & lf2r_state_seq{ns}(2:end) == 2);
            cur_down_trans = new_seg_inds(ns,1)-1+find(lf2r_state_seq{ns}(1:end-1) == 2 & lf2r_state_seq{ns}(2:end) == 1);
            cur_up_trans(cur_up_trans > cur_down_trans(end)) = [];
            cur_down_trans(cur_down_trans < cur_up_trans(1)) = [];
            
            n_lf2r_ups = length(cur_up_trans);
            
            %get rid of LF2 up transitions where LF8 was desynched
            bad_states = find(is_desynch8(cur_up_trans) == 1);
            cur_up_trans(bad_states) = [];
            cur_down_trans(bad_states) = [];
            
            for i = 1:length(cur_up_trans)
                
                %find nearest lfp up transition
                [dummy,near_lfp_up] = min(abs(cur_up_trans(i) - ov_up_trans8));
                
                if n_lfp_ups > near_lfp_up+1
                    %lf2r and nearest lfp up-transition indices
                    lf2r_up_inds{d} = [lf2r_up_inds{d} cur_up_trans(i)];
                    lf2r_down_inds{d} = [lf2r_down_inds{d} cur_down_trans(i)];
                    corresp_lfp_ind{d} = [corresp_lfp_ind{d} near_lfp_up];
                    
                    %lf2r up state dur
                    lf2r_updur{d} = [lf2r_updur{d} (cur_down_trans(i)-cur_up_trans(i))/Fsd];
                    
                    %store the indices of any skipped lfp up transitions
                    skipped_lfp_ups = find(ov_up_trans8 > ov_up_trans8(near_lfp_up) ...
                        & ov_up_trans8 < cur_down_trans(i));
                    lf2r_skipped_lfpupind{d} = [lf2r_skipped_lfpupind{d}; skipped_lfp_ups];
                    
                    %store the LFP up-transition lag in sec
                    cur_lf2r_up_lag = (cur_up_trans(i)-ov_up_trans8(near_lfp_up));
                    lf2r_up_lag{d} = [lf2r_up_lag{d} cur_lf2r_up_lag/Fsd];
                    
                    %store lf2r up state dur in lfp cycles
                    cur_lfpc = ov_lf8_period_f(cur_down_trans(i) - cur_lf2r_up_lag) ...
                        - ov_lf8_period_f(cur_up_trans(i)-cur_lf2r_up_lag);
                    lf2r_updur_lfpc{d} = [lf2r_updur_lfpc{d} cur_lfpc];
                    
                    %find corresponding down transition
                    next_lfp_down = ov_down_trans8(near_lfp_up);
                    if cur_down_trans(i) < next_lfp_down
                        corresp_lfp_down = next_lfp_down;
                    else
                        corresp_lfp_down = ov_down_trans8(find(ov_down_trans8 < cur_down_trans(i),1,'last'));
                    end
                    if ~isempty(corresp_lfp_down)
                        cur_lf2r_down_lag = cur_down_trans(i)-corresp_lfp_down;
                        lf2r_down_lag{d} = [lf2r_down_lag{d} cur_lf2r_down_lag/Fsd];
                    else
                        lf2r_down_lag{d} = [lf2r_down_lag{d} nan];
                    end
                    
                    corresp_lfp_updur{d} = [corresp_lfp_updur{d} (ov_down_trans8(near_lfp_up)-ov_up_trans8(near_lfp_up))/Fsd];
                    corresp_lfp_downdur{d} = [corresp_lfp_downdur{d} (ov_up_trans8(near_lfp_up+1)-ov_down_trans8(near_lfp_up))/Fsd];
                    corresp_lfp_cycledur{d} = [corresp_lfp_cycledur{d} (ov_up_trans8(near_lfp_up+1)-ov_up_trans8(near_lfp_up))/Fsd];
                    corresp_lfp_dutycycle{d} = [corresp_lfp_dutycycle{d} corresp_lfp_updur{d}(end)/corresp_lfp_cycledur{d}(end)];
                    %                 corresp_lfp_amp{d} = [corresp_lfp_amp{d} lf8_uds_amp{ns}(near_lfp_up)];
                    
                end
            end
            
            %             for i = 1:length(cur_up_trans8)
            %                 %find nearest lf2r down transition
            %                 [dummy,near_lf2r_down] = min(abs(cur_up_trans8(i) - cur_down_trans));
            %                 lup_mdown_lag{d} = [lup_mdown_lag{d}; cur_up_trans8(i)-cur_down_trans(near_lf2r_down)];
            %             end
            
        end
        
        pers_fract_within_cycle_lf2r(d) = length(find(lf2r_updur_lfpc{d}>lf8_up_fract))/length(lf2r_updur_lfpc{d});
        pers_fract_across_cycles_lf2r(d) = length(find(lf2r_updur_lfpc{d} > 1))/length(lf2r_updur_lfpc{d});
        pers_fract_within_np_lf2r(d) = length(find(lf2r_updur_lfpc{d} >lf8_up_fract & lf2r_updur_lfpc{d} < 1))...
            /length(find(lf2r_updur_lfpc{d} < 1));
        persistent_dur_lf2r(d) = nanmean(lf2r_updur{d}(lf2r_updur_lfpc{d} > 1));
        mean_up_lag_lf2r(d) = nanmean(lf2r_up_lag{d});
        mean_down_lag_lf2r(d) = nanmean(lf2r_down_lag{d});
        median_up_lag_lf2r(d) = nanmedian(lf2r_up_lag{d});
        median_down_lag_lf2r(d) = nanmedian(lf2r_down_lag{d});
        
        clear lf8_period*
    end
    
    %%
%     if bad_lf3(d) == 0 && uds_time_3 > 250
%         
%         datalen = length(lf3);
%         load ec_hmm_state_seq3
%         lf3_state_seq = hmm_bbstate_seq3;
%         
%         [new_seg_inds] = resample_uds_seg_inds(hmm3.UDS_segs,hmm3.Fs,Fsd,length(lf3));
%         seg_durs = new_seg_inds(:,2)-new_seg_inds(:,1)+1;
%         
%         lf3_up_inds{d} = [];
%         lf3_down_inds{d} = [];
%         lf3_updur{d} = [];
%         lfp_up_lag{d} = [];
%         lfp_down_lag{d} = [];
%         lf3_updur_lfpc{d} = [];
%         corresp_lfp_updur{d} = [];
%         corresp_lfp_downdur{d} = [];
%         corresp_lfp_cycledur{d} = [];
%         corresp_lfp_dutycycle{d} = [];
%         corresp_lfp_alf3{d} = [];
%         corresp_lfp_ind{d} = [];
%         lf3_skipped_lfpupind{d} = [];
%         lup_mdown_lag{d} = [];
%         
%         for ns = 1:hmm3.Nsegs
%             
%             cur_up_trans = new_seg_inds(ns,1)-1+find(lf3_state_seq{ns}(1:end-1) == 1 & lf3_state_seq{ns}(2:end) == 2);
%             cur_down_trans = new_seg_inds(ns,1)-1+find(lf3_state_seq{ns}(1:end-1) == 2 & lf3_state_seq{ns}(2:end) == 1);
%             cur_up_trans(cur_up_trans > cur_down_trans(end)) = [];
%             cur_down_trans(cur_down_trans < cur_up_trans(1)) = [];
%             
%             n_lf3_ups = length(cur_up_trans);
%             
%             %get rid of LF2 up transitions where LF8 was desynched
%             bad_states = find(is_desynch8(cur_up_trans) == 1);
%             cur_up_trans(bad_states) = [];
%             cur_down_trans(bad_states) = [];
%             
%             for i = 1:length(cur_up_trans)
%                 
%                 %find nearest lfp up transition
%                 [dummy,near_lfp_up] = min(abs(cur_up_trans(i) - ov_up_trans8));
%                 
%                 if n_lfp_ups > near_lfp_up+1
%                     %lf3 and nearest lfp up-transition indices
%                     lf3_up_inds{d} = [lf3_up_inds{d} cur_up_trans(i)];
%                     lf3_down_inds{d} = [lf3_down_inds{d} cur_down_trans(i)];
%                     corresp_lfp_ind{d} = [corresp_lfp_ind{d} near_lfp_up];
%                     
%                     %lf3 up state dur
%                     lf3_updur{d} = [lf3_updur{d} (cur_down_trans(i)-cur_up_trans(i))/Fsd];
%                     
%                     %store the indices of any skipped lfp up transitions
%                     skipped_lfp_ups = find(ov_up_trans8 > ov_up_trans8(near_lfp_up) ...
%                         & ov_up_trans8 < cur_down_trans(i));
%                     lf3_skipped_lfpupind{d} = [lf3_skipped_lfpupind{d}; skipped_lfp_ups];
%                     
%                     %store the LFP up-transition lag in sec
%                     cur_lfp_up_lag = (cur_up_trans(i)-ov_up_trans8(near_lfp_up));
%                     lfp_up_lag{d} = [lfp_up_lag{d} cur_lfp_up_lag/Fsd];
%                     
%                     %store lf3 up state dur in lfp cycles
%                     cur_lfpc = ov_lf8_period_f(cur_down_trans(i) - cur_lfp_up_lag) ...
%                         - ov_lf8_period_f(cur_up_trans(i)-cur_lfp_up_lag);
%                     lf3_updur_lfpc{d} = [lf3_updur_lfpc{d} cur_lfpc];
%                     
%                     %find corresponding down transition
%                     next_lfp_down = ov_down_trans8(near_lfp_up);
%                     if cur_down_trans(i) < next_lfp_down
%                         corresp_lfp_down = next_lfp_down;
%                     else
%                         corresp_lfp_down = ov_down_trans8(find(ov_down_trans8 < cur_down_trans(i),1,'last'));
%                     end
%                     if ~isempty(corresp_lfp_down)
%                         cur_lfp_down_lag = cur_down_trans(i)-corresp_lfp_down;
%                         lfp_down_lag{d} = [lfp_down_lag{d} cur_lfp_down_lag/Fsd];
%                     else
%                         lfp_down_lag{d} = [lfp_down_lag{d} nan];
%                     end
%                     
%                     corresp_lfp_updur{d} = [corresp_lfp_updur{d} (ov_down_trans8(near_lfp_up)-ov_up_trans8(near_lfp_up))/Fsd];
%                     corresp_lfp_downdur{d} = [corresp_lfp_downdur{d} (ov_up_trans8(near_lfp_up+1)-ov_down_trans8(near_lfp_up))/Fsd];
%                     corresp_lfp_cycledur{d} = [corresp_lfp_cycledur{d} (ov_up_trans8(near_lfp_up+1)-ov_up_trans8(near_lfp_up))/Fsd];
%                     corresp_lfp_dutycycle{d} = [corresp_lfp_dutycycle{d} corresp_lfp_updur{d}(end)/corresp_lfp_cycledur{d}(end)];
%                     %                 corresp_lfp_amp{d} = [corresp_lfp_amp{d} lf8_uds_amp{ns}(near_lfp_up)];
%                     
%                 end
%             end
%             
%             %             for i = 1:length(cur_up_trans8)
%             %                 %find nearest lf3 down transition
%             %                 [dummy,near_lf3_down] = min(abs(cur_up_trans8(i) - cur_down_trans));
%             %                 lup_mdown_lag{d} = [lup_mdown_lag{d}; cur_up_trans8(i)-cur_down_trans(near_lf3_down)];
%             %             end
%             
%         end
%         
%         pers_fract_within_cycle_lf3(d) = length(find(lf3_updur_lfpc{d}>lf8_up_fract))/length(lf3_updur_lfpc{d});
%         pers_fract_across_cycles_lf3(d) = length(find(lf3_updur_lfpc{d} > 1))/length(lf3_updur_lfpc{d});
%         pers_fract_within_np_lf3(d) = length(find(lf3_updur_lfpc{d} >lf8_up_fract & lf3_updur_lfpc{d} < 1))...
%             /length(find(lf3_updur_lfpc{d} < 1));
%         persistent_dur_lf3(d) = nanmean(lf3_updur{d}(lf3_updur_lfpc{d} > 1));
%         mean_up_lag_lf3(d) = nanmean(lfp_up_lag{d});
%         mean_down_lag_lf3(d) = nanmean(lfp_down_lag{d});
%         median_up_lag_lf3(d) = nanmedian(lfp_up_lag{d});
%         median_down_lag_lf3(d) = nanmedian(lfp_down_lag{d});
%         
%         clear lf8_period*
%     end
end

cd F:\WC_Germany\overall_EC\
save corresponding_lfp_state_data_hipp_lfp pers_fract* mean*lag* median*lag* lf*lag* lf*_up* lf*_down* corresp* lf*_skipped_lfpupind lup*

