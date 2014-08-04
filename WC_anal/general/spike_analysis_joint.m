clear all
close all

load C:\WC_Germany\Persistent_activity\dir_tree_update
load C:\WC_Germany\Persistent_activity\UDS_dur_raw\UDS_raw_data
load C:\WC_Germany\Persistent_activity\UDS_dist\synch_UDS_inds
load C:\WC_Germany\Persistent_activity\lf8_period_f_data2

Fs = 2016;
dsf = 8;
Fsd = Fs/dsf;
forwardlag = 5*Fsd;
backlag = 1*Fsd;
binsize = 50;
lags = -backlag:binsize:forwardlag;
phase_hist = linspace(0,360,15);

for d = 1:length(dir_array)

    cd(dir_array{d})
    disp(['session ' num2str(d)])

    load spike_time_jmm
    load used_data wcv lf8
    wcv_d = downsample(wcv,dsf);
    datalen = length(wcv_d);
    spkids = round(spkid/dsf);

    tsu_spike_mat = zeros(length(synch_ups{d}),length(lags));
    tsu_phase_mat = zeros(length(synch_ups{d}),length(lags));
    joint_spike_mat = zeros(length(phase_hist),length(lags));
    tsu_time_mat = zeros(length(synch_ups{d}),length(lags));
    joint_time_mat = zeros(length(phase_hist),length(lags));
    
    for i = 1:length(synch_ups{d})
        cur_up_trans = up_trans{d}(synch_ups{d}(i));
        cur_down_trans = down_trans{d}(synch_ups{d}(i));
        if cur_up_trans > backlag && (datalen-cur_up_trans) > forwardlag
            cur_spikes = spkids(find(spkids > (cur_up_trans - backlag) & spkids < (cur_up_trans + forwardlag)));
            spike_cnt = hist(cur_spikes-cur_up_trans,lags);
            up_trans_mat(i,:) = spike_cnt;
            tsu_time_mat(i,:) = ones(size(lags));
            %don't count times when mp was in down state
            if cur_down_trans - cur_up_trans < forwardlag
                state_dur = cur_down_trans-cur_up_trans;
                tsu_time_mat(i,backlag + state_dur:end) = 0;
            end
            %now extract phase at each time
            if length(lf8_period_p{d}) > cur_up_trans+forwardlag
                cur_phases = lf8_period_p{d}(cur_up_trans-backlag:cur_up_trans+forwardlag);
                lf8_phase_mat(i,:) = cur_phases;
                %don't count tmies when mp was in down state
                if cur_down_trans - cur_up_trans < forwardlag
                     state_dur = cur_down_trans-cur_up_trans;
                     lf8_phase_mat(i,backlag+state_dur:end) = nan;
                end
                
                % cycle through all spikes
                rel_spike_times = cur_spikes - cur_up_trans;
                spike_phases = lf8_period_p{d}(cur_spikes);
                for s = 1:length(cur_spikes)
                    time_bin = find(lags > rel_spike_times(s),1,'first');
                    phase_bin = find(phase_hist > spike_phases(s),1,'first');
                    joint_spike_mat(phase_bin,time_bin) = joint_spike_mat(phase_bin,time_bin)+1;
                end
                
                %cycle through all lags and find time spent at each phase
                %range
                for l = 1:(length(lags)-1)
                    beg_id = (l-1)*binsize+1;
                    end_id = l*binsize;
                    phase_binned = hist(cur_phases(beg_id:end_id),phase_hist);
                    phase_binned = phase_binned';
                    joint_time_mat(:,l) = joint_time_mat(:,l) + phase_binned;
                end
                    
                
            else
                lf8_phase_mat(i,:) = nan;
            end
        else
            up_trans_mat(i,:) = nan;
            lf8_phase_mat(i,:) = nan;
        end
    end

    joint_rate = joint_spike_mat./joint_time_mat;
    joint_rate(isnan(joint_rate)) = 0;
    joint_rate(isinf(joint_rate)) = 0;
    joint_rate = joint_rate*Fsd;
    sm_joint_rate = gauss_smooth2(joint_rate,2,9);
    pcolor(lags/Fsd,phase_hist,joint_rate);shading flat
    caxis([0 20]);colorbar
    t_names = ['C:\WC_Germany\Persistent_activity\spike_analysis\joint_rate_' f_names{d}];
        xlabel('Time since up trans (s)')
        ylabel('Phase (degrees)')
        colorbar
        print('-dpng',t_names)
        close

end

