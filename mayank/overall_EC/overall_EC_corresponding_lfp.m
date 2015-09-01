clear all

load F:\WC_Germany\overall_EC\overall_EC_dir.mat

drive_letter = 'F';
dsf = 8;
Fsd = 2016/dsf;
lf8_up_fract = 0.5;

for d = 1:length(sess_data)
    
    fprintf('session %d\n',d)
    cdir = sess_data(d).directory;
    cdir(1) = drive_letter;
    cd(cdir)
    
    load ./used_data wcv
    wcv_d = downsample(wcv,dsf);
    datalen = length(wcv_d);
    load ./ec_hmm_state_seq
    mp_state_seq = hmm_bbstate_seq;
    load ./ec_hmm_state_seq8
    lfp_state_seq = hmm_bbstate_seq8;
   
    load ./lf8_period_data
    
    [new_seg_inds] = resample_uds_seg_inds(hmm.UDS_segs,hmm.Fs,Fsd,length(wcv_d));
    seg_durs = new_seg_inds(:,2)-new_seg_inds(:,1)+1;
    
    mp_up_inds{d} = [];
    mp_down_inds{d} = [];
    mp_updur{d} = [];
    lfp_up_lag{d} = [];
    lfp_down_lag{d} = [];
    mp_updur_lfpc{d} = [];
    corresp_lfp_updur{d} = [];
    corresp_lfp_downdur{d} = [];
    corresp_lfp_cycledur{d} = [];
    corresp_lfp_dutycycle{d} = [];
    corresp_lfp_amp{d} = [];
    corresp_lfp_ind{d} = [];
    mp_skipped_lfpupind{d} = [];
    lup_mdown_lag{d} = [];

    for ns = 1:hmm.Nsegs
        
        cur_up_trans = find(mp_state_seq{ns}(1:end-1) == 1 & mp_state_seq{ns}(2:end) == 2);
        cur_down_trans = find(mp_state_seq{ns}(1:end-1) == 2 & mp_state_seq{ns}(2:end) == 1);
        cur_up_trans(cur_up_trans > cur_down_trans(end)) = [];
        cur_down_trans(cur_down_trans < cur_up_trans(1)) = [];
%         cur_up_trans = cur_up_trans + new_seg_inds(ns,1)-1;
%         cur_down_trans = cur_down_trans + new_seg_inds(ns,1)-1;       

        cur_up_trans8 = find(lfp_state_seq{ns}(1:end-1) == 1 & lfp_state_seq{ns}(2:end) == 2);
        cur_down_trans8 = find(lfp_state_seq{ns}(1:end-1) == 2 & lfp_state_seq{ns}(2:end) == 1);
        cur_up_trans8(cur_up_trans8 > cur_down_trans8(end)) = [];
        cur_down_trans8(cur_down_trans8 < cur_up_trans8(1)) = [];
%         cur_up_trans8 = cur_up_trans8 + new_seg_inds(ns,1)-1;
%         cur_down_trans8 = cur_down_trans8 + new_seg_inds(ns,1)-1;
        
        n_lfp_ups = length(cur_up_trans8);
        n_mp_ups = length(cur_up_trans);
        
        for i = 1:length(cur_up_trans)
            
            %find nearest lfp up transition
            [dummy,near_lfp_up] = min(abs(cur_up_trans(i) - cur_up_trans8));
            
            if n_lfp_ups > near_lfp_up+1
                %mp and nearest lfp up-transition indices
                mp_up_inds{d} = [mp_up_inds{d} (new_seg_inds(ns,1)-1+cur_up_trans(i))];
                mp_down_inds{d} = [mp_down_inds{d} (new_seg_inds(ns,1)-1+cur_down_trans(i))];
                corresp_lfp_ind{d} = [corresp_lfp_ind{d} near_lfp_up];
                
                %mp up state dur
                mp_updur{d} = [mp_updur{d} (cur_down_trans(i)-cur_up_trans(i))/Fsd];
                
                %store the indices of any skipped lfp up transitions
                skipped_lfp_ups = find(cur_up_trans8 > cur_up_trans8(near_lfp_up) ...
                    & cur_up_trans8 < cur_down_trans(i));
                mp_skipped_lfpupind{d} = [mp_skipped_lfpupind{d}; skipped_lfp_ups];
                
                %store the LFP up-transition lag in sec
                cur_lfp_up_lag = (cur_up_trans(i)-cur_up_trans8(near_lfp_up));
                lfp_up_lag{d} = [lfp_up_lag{d} cur_lfp_up_lag/Fsd];
                
                %store mp up state dur in lfp cycles
                cur_lfpc = lf8_period_f{ns}(cur_down_trans(i) - cur_lfp_up_lag) ...
                    - lf8_period_f{ns}(cur_up_trans(i)-cur_lfp_up_lag);
                mp_updur_lfpc{d} = [mp_updur_lfpc{d} cur_lfpc];
                
                %find corresponding down transition
                next_lfp_down = cur_down_trans8(near_lfp_up);
                if cur_down_trans(i) < next_lfp_down
                    corresp_lfp_down = next_lfp_down;
                else
                    corresp_lfp_down = cur_down_trans8(find(cur_down_trans8 < cur_down_trans(i),1,'last'));
                end
                if ~isempty(corresp_lfp_down)
                    cur_lfp_down_lag = cur_down_trans(i)-corresp_lfp_down;
                    lfp_down_lag{d} = [lfp_down_lag{d} cur_lfp_down_lag/Fsd];
                else
                    lfp_down_lag{d} = [lfp_down_lag{d} nan];
                end
                
                corresp_lfp_updur{d} = [corresp_lfp_updur{d} (cur_down_trans8(near_lfp_up)-cur_up_trans8(near_lfp_up))/Fsd];
                corresp_lfp_downdur{d} = [corresp_lfp_downdur{d} (cur_up_trans8(near_lfp_up+1)-cur_down_trans8(near_lfp_up))/Fsd];
                corresp_lfp_cycledur{d} = [corresp_lfp_cycledur{d} (cur_up_trans8(near_lfp_up+1)-cur_up_trans8(near_lfp_up))/Fsd];
                corresp_lfp_dutycycle{d} = [corresp_lfp_dutycycle{d} corresp_lfp_updur{d}(end)/corresp_lfp_cycledur{d}(end)];
%                 corresp_lfp_amp{d} = [corresp_lfp_amp{d} lf8_uds_amp{ns}(near_lfp_up)];
                
            end
        end
        
        for i = 1:length(cur_up_trans8)
            %find nearest MP down transition
            [dummy,near_mp_down] = min(abs(cur_up_trans8(i) - cur_down_trans));
            lup_mdown_lag{d} = [lup_mdown_lag{d}; cur_up_trans8(i)-cur_down_trans(near_mp_down)];
        end
        
    end
        
    pers_fract_within_cycle(d) = length(find(mp_updur_lfpc{d}>lf8_up_fract))/length(mp_updur_lfpc{d});
    pers_fract_across_cycles(d) = length(find(mp_updur_lfpc{d} > 1))/length(mp_updur_lfpc{d});
    pers_fract_within_np(d) = length(find(mp_updur_lfpc{d} >lf8_up_fract & mp_updur_lfpc{d} < 1))...
        /length(find(mp_updur_lfpc{d} < 1));
    persistent_dur(d) = nanmean(mp_updur{d}(mp_updur_lfpc{d} > 1));
    mean_up_lag(d) = nanmean(lfp_up_lag{d});
    mean_down_lag(d) = nanmean(lfp_down_lag{d});
    median_up_lag(d) = nanmedian(lfp_up_lag{d});
    median_down_lag(d) = nanmedian(lfp_down_lag{d});
    
    clear lf8_period*
end


cd F:\WC_Germany\overall_EC\
save corresponding_lfp_state_data pers_fract* mean*lag median*lag lfp*lag mp_up* mp_down* corresp* mp_skipped_lfpupind lup*

