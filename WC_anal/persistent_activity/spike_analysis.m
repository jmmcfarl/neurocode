clear all
close all

load C:\WC_Germany\Persistent_activity\dir_tree_update
load C:\WC_Germany\Persistent_activity\UDS_dur_raw\UDS_raw_data
load C:\WC_Germany\Persistent_activity\UDS_dist\synch_UDS_inds
load C:\WC_Germany\Persistent_activity\lf8_period_f_data2

Fs = 2016;
dsf = 8;
Fsd = Fs/dsf;
forwardlag = 10*Fsd;
backlag = 1*Fsd;
binsize = 20;
lags = -backlag:binsize:forwardlag;
phase_hist = linspace(0,360,30);

for d = 1:length(dir_array)

    cd(dir_array{d})
    disp(['session ' num2str(d)])

    load spike_time_jmm
    load used_data wcv lf8
    wcv_d = downsample(wcv,dsf);
    datalen = length(wcv_d);
    %    plot(wcv)
    %    hold on
    %    plot(spkid,wcv(spkid),'r.')

    spkids = round(spkid/dsf);

    %     wcv_isup = zeros(size(wcv_d));
    %     for i = 1:length(synch_ups{d})
    %        wcv_isup(up_trans{d}(synch_ups{d}(i)):down_trans{d}(synch_ups{d}(i))) = 1;
    %     end

    up_trans_mat = zeros(length(synch_ups{d}),length(lags));
    lf8_phase_mat = zeros(length(synch_ups{d}),length(lags));
    for i = 1:length(synch_ups{d})
        cur_up_trans = up_trans{d}(synch_ups{d}(i));
        cur_down_trans = down_trans{d}(synch_ups{d}(i));
        if cur_up_trans > backlag && (datalen-cur_up_trans) > forwardlag
            cur_spikes = spkids(find(spkids > (cur_up_trans - backlag) & spkids < (cur_up_trans + forwardlag)));
            up_trans_mat(i,:) = hist(cur_spikes - cur_up_trans,lags);
            if cur_down_trans - cur_up_trans < forwardlag
                state_dur = cur_down_trans-cur_up_trans;
                up_trans_mat(i,backlag + state_dur:end) = nan;
            end
            if length(lf8_period_p{d}) > cur_up_trans+forwardlag
                lf8_phase_mat(i,:) = lf8_period_p{d}(cur_up_trans-backlag:binsize:cur_up_trans+forwardlag);
                if cur_down_trans - cur_up_trans < forwardlag
                     state_dur = cur_down_trans-cur_up_trans;
                     lf8_phase_mat(i,backlag+state_dur:end) = nan;
                end
            else
                lf8_phase_mat(i,:) = nan;
            end
        else
            up_trans_mat(i,:) = nan;
            lf8_phase_mat(i,:) = nan;
        end
    end

%       %% compute spike probability contingent on an LFP up transition that
%     %% occurs during an MP up state  
%     up_cnts = 1;
%     clear cond_up_trans_mat
%     for i = 1:length(synch_ups{d})
%         cur_up_trans = up_trans{d}(synch_ups{d}(i));
%         cur_down_trans = down_trans{d}(synch_ups{d}(i));
%         %find all lfp up transitions occuring during MP up state and at least 100ms after MP up transition
%         good_lfp_ups = up_trans8{d}(find(up_trans8{d} > (cur_up_trans + round(0.1*Fsd)) & up_trans8{d} < (cur_down_trans - round(0.1*Fsd))));
%         for cu = 1:length(good_lfp_ups)
%             next_lfp_up = good_lfp_ups(cu);
%             if next_lfp_up > backlag & (datalen-next_lfp_up) > forwardlag
%                 cur_spikes = spkids(find(spkids > (next_lfp_up - backlag) & spkids < (next_lfp_up + forwardlag)));
%                 cond_up_trans_mat(up_cnts,:) = hist(cur_spikes - next_lfp_up,lags);
%                 %set values before mp up transition to nans
%                 lfp_delay = next_lfp_up-cur_up_trans;
%                 if lfp_delay < backlag
%                    cond_up_trans_mat(up_cnts,1:backlag-lfp_delay) = nan; 
%                 end
%                 %set values after mp down transition to nans
%                 mp_delay = cur_down_trans - next_lfp_up;
%                 if mp_delay < forwardlag
%                    cond_up_trans_mat(up_cnts,backlag+length(lags)-mp_delay:end) = nan; 
%                 end
%                 up_cnts = up_cnts+1;
%             end
%         end
%     end

    lf8_cond_spkmat = zeros(length(synch_ups{d}),length(phase_hist));
    time_sls = [];
    lf8_phase = [];
    spike_vec = [];
    for i = 1:length(synch_ups{d})
        cur_up_trans = up_trans{d}(synch_ups{d}(i));
        cur_down_trans = down_trans{d}(synch_ups{d}(i));
        if cur_down_trans-cur_up_trans > Fsd
            cur_up_trans = cur_up_trans+Fsd;
            
        if length(lf8_period_p{d}) > cur_down_trans
            cur_lf8_phase = lf8_period_p{d}(cur_up_trans:cur_down_trans);
            cur_spk_ids = spkids(find(spkids > cur_up_trans & spkids < cur_down_trans)) - cur_up_trans;
           
            spike_phase_cnt = hist(cur_lf8_phase(cur_spk_ids),phase_hist);
            total_phase_cnt = hist(cur_lf8_phase,phase_hist);
            total_phase_cnt(total_phase_cnt == 0) = nan;
            lf8_cond_spkmat(i,:) = spike_phase_cnt./total_phase_cnt;
            time_sls = [time_sls (0:(cur_down_trans-cur_up_trans))/Fsd];
            lf8_phase = [lf8_phase lf8_period_p{d}(cur_up_trans:cur_down_trans)];
                cur_spike_vec = zeros(1,cur_down_trans-cur_up_trans+1);
                cur_spike_vec(cur_spk_ids) = 1;
                spike_vec = [spike_vec cur_spike_vec];
            
            
        else
            lf8_cond_spkmat(i,:) = nan;
        end
        else
            lf8_cond_spkmat(i,:) = nan;
        end
    end

    response = spike_vec;
    predictor = [sin(lf8_phase*pi/180); time_sls];
    
    num_good_trials = sum(~isnan(up_trans_mat));
%     num_good_trials_ph = sum(~isnan(lf8_phase_mat));
    for i = 1:length(lags)
        if num_good_trials(i) < 10
            up_trans_mat(:,i) = nan;
        end
%         if num_good_trials_ph(i) < 10
%             lf8_phase_mat(:,i) = nan;
%         end
    end
    
    %convert to units of rate
    up_trans_mat = up_trans_mat*Fsd/binsize;
    lf8_cond_spkmat = lf8_cond_spkmat*Fsd;
%     cond_up_trans_mat = cond_up_trans_mat*Fsd/binsize;
    
    avg_rate_traj(d,:) = jmm_smooth_1d(nanmean(up_trans_mat),4);
    avg_phase_traj(d,:) = jmm_smooth_1d(nanmean(lf8_phase_mat),4);
    avg_lfp_response(d,:) = jmm_smooth_1d(nanmean(cond_up_trans_mat),4);
    avg_lf8_cond_traj(d,:) = nanmean(lf8_cond_spkmat);
    num_used_trials = sum(~isnan(lf8_cond_spkmat));
    se_lf8_cond_traj(d,:) = nanstd(lf8_cond_spkmat)./sqrt(num_used_trials);

%         errorbar(phase_hist,avg_lf8_cond_traj(d,:),se_lf8_cond_traj(d,:))
%          t_names = ['C:\WC_Germany\Persistent_activity\spike_analysis\lf8_cond_spike_' f_names{d}];
%        xlabel('LF8 UDS phase (degrees)')
%        ylabel('Average firing rate(Hz)')
%        xlim([0 360])
%        print('-dpng',t_names)
%        close
%     
%         t_names = ['C:\WC_Germany\Persistent_activity\spike_analysis\up_trans_spike_' f_names{d}];
%         plot(lags/Fsd,avg_rate_traj(d,:))
%         xlabel('Time since up trans (s)')
%         ylabel('Average rate (Hz)')
%     
%         print('-dpng',t_names)
%         close


end

%compute overall averages
for i = 1:17
    avg_rate_traj(i,:) = avg_rate_traj(i,:)/max(avg_rate_traj(i,:));
    avg_lf8_cond_traj(i,:) = avg_lf8_cond_traj(i,:)/nanmax(avg_lf8_cond_traj(i,:));
    avg_lfp_response(i,:) = avg_lfp_response(i,:)/nanmax(avg_lfp_response(i,:));
end

% figure(1)
% errorbar(lags/Fsd,nanmean(avg_rate_traj),nanstd(avg_rate_traj)/sqrt(17))
% 
% figure(2)
% errorbar(phase_hist,nanmean(avg_lf8_cond_traj),nanstd(avg_lf8_cond_traj)/sqrt(17))

