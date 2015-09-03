clear all
% close all

load C:\WC_Germany\Persistent_activity\dir_tree_update
load C:\WC_Germany\Persistent_activity\UDS_dur_raw\UDS_raw_data
load C:\WC_Germany\Persistent_activity\UDS_dist\synch_UDS_inds
load C:\WC_Germany\Persistent_activity\lf8_period_f_data2

Fs = 2016;
dsf = 8;
Fsd = Fs/dsf;
forwardlag = 5*Fsd;
backlag = 0;
binsize = 1;
lags = -backlag:binsize:forwardlag;
phase_hist = linspace(0,360,30);

for d = 1:length(dir_array)
    cd(dir_array{d})
    disp(['session ' num2str(d)])

    load spike_time_jmm
    load used_data wcv lf8
    wcv_d = downsample(wcv,dsf);
    datalen = length(wcv_d);

    spkids = round(spkid/dsf);

    
    net_spike_phase = zeros(1,length(phase_hist));
    net_phase = zeros(1,length(phase_hist));
    net_spike_tslu = zeros(1,length(lags));
    net_tslu = zeros(1,length(lags));
    time_sls = [];
    lf8_phase = [];
    spike_vec = [];
    lf8_phase_mat = nan(length(synch_ups{d}),length(lags));
    cur_spike_vec = zeros(1,length(lags));

    for i = 1:length(synch_ups{d})
        
        cur_up_trans = up_trans{d}(synch_ups{d}(i))+round(0.1*Fsd);
        cur_down_trans = down_trans{d}(synch_ups{d}(i))-round(0.1*Fsd);
        
        if length(lf8_period_p{d}) > cur_down_trans
            cur_spk_ids = spkids(find(spkids > cur_up_trans & spkids < cur_down_trans)) - cur_up_trans+1;
            cur_spk_ids(cur_spk_ids >=forwardlag) = []; %don't count spikes after a max lag
            net_spike_tslu = net_spike_tslu+hist(cur_spk_ids,lags); %increment net spike count
            
            %increment count of time spent at each tslu
            if cur_down_trans-cur_up_trans > forwardlag
                net_tslu = net_tslu + 1;
            else
                end_lag_bin = round((cur_down_trans-cur_up_trans)/binsize)+1;
                net_tslu(1:end_lag_bin) = net_tslu(1:end_lag_bin) + 1;
            end
            
            cur_lf8_phase = lf8_period_p{d}(cur_up_trans:cur_down_trans);
            
            %don't take phases beyond max lag
            if length(cur_lf8_phase) > forwardlag
                cur_lf8_phase(forwardlag+2:end) = [];
                lf8_phase_mat(i,:) = cur_lf8_phase;
            else
                lf8_phase_mat(i,1:length(cur_lf8_phase)) = cur_lf8_phase;
            end
            
            spike_phase_cnt = hist(cur_lf8_phase(cur_spk_ids),phase_hist);
            total_phase_cnt = hist(cur_lf8_phase,phase_hist);
            
            net_spike_phase = net_spike_phase + spike_phase_cnt;
            net_phase = net_phase + total_phase_cnt;

            if cur_down_trans-cur_up_trans > forwardlag
                time_sls = [time_sls 0:forwardlag];
            else
                time_sls = [time_sls 0:(cur_down_trans-cur_up_trans)];
            end
            lf8_phase = [lf8_phase cur_lf8_phase];
            
            
            cur_spike_vec(cur_spk_ids) = 1;
            spike_vec = [spike_vec cur_spike_vec];
               
        else
            lf8_phase_mat(i,:) = nan;
        end
    end
    
    %calculate mean rate
    mean_rate(d) = sum(net_spike_tslu)/sum(net_tslu)*Fsd;
    
    %calculate empirical quantities
    emp_phase_rate(d,:) = net_spike_phase./net_phase*Fsd;
    
    %smooth empirical phase rate slightly
%     emp_phase_rate(d,:) = jmm_smooth_1d(emp_phase_rate(d,:),2);
    
    emp_tslu_rate(d,:) = jmm_smooth_1d(net_spike_tslu./net_tslu*Fsd,20);    
    for c = 1:length(lags)
       used_phases = lf8_phase_mat(:,c);
       used_phases(isnan(used_phases)) = [];
       if ~isempty(used_phases)
            emp_phase_tslu(d,c) = circ_mean(used_phases*pi/180);
       else
           emp_phase_tslu(d,c) = nan;
       end
    end
    emp_phase_tslu(d,:) = emp_phase_tslu(d,:)*180/pi;
    neg_vals = find(emp_phase_tslu(d,:) < 0);
    emp_phase_tslu(d,neg_vals) = emp_phase_tslu(d,neg_vals) + 360;
    
    
    [boot_phase_rate(d,:),phase_time_dist] = get_boot_phase_rate_v2(lf8_phase,time_sls,phase_hist,lags,emp_tslu_rate(d,:),500);
   
    
    clear net_spike_phase net_pahse lf8_phase time_sls 

%     figure(1)
%     plot(lags/Fsd,emp_tslu_rate(d,:))
%     t_names = ['C:\WC_Germany\Persistent_activity\spike_analysis\rate_tsu_' f_names{d}];
%     xlabel('Time since up trans (s)')
%     ylabel('Mean Rate (Hz)')
%     print('-dpng',t_names)
%     close
% 
%     figure(2)
%     plot(lags/Fsd,emp_phase_tslu(d,:))
%     t_names = ['C:\WC_Germany\Persistent_activity\spike_analysis\phase_tsu_' f_names{d}];
%     xlabel('Time since up trans (s)')
%     ylabel('Circular Mean Phase (degrees)')
%     print('-dpng',t_names)
%     close
% 
%     figure(3)
%     plot(phase_hist,boot_phase_rate(d,:))
%     hold on
%     plot(phase_hist,emp_phase_rate(d,:),'r')
% %     errorbar(phase_hist,emp_phase_rate(d,:),emp_phase_rate_unc(d,:),'r')
%     legend('Independent','Empirical')
%     t_names = ['C:\WC_Germany\Persistent_activity\spike_analysis\phase_rate_' f_names{d}];
%     xlabel('Phase (degrees)')
%     ylabel('Mean Rate (Hz)')
%     print('-dpng',t_names)
%     close
%     
%     figure(4)
%     pcolor(lags/Fsd,phase_hist,phase_time_dist)   
%     shading flat
%     t_names = ['C:\WC_Germany\Persistent_activity\spike_analysis\phase_time_dist_' f_names{d}];
%     xlabel('Time since up trans (s)')
%     ylabel('Circular Mean Phase (degrees)')
%     print('-dpng',t_names)
%     close

    
%     %normalize to max rate
%     emp_phase_rate(d,:) = emp_phase_rate(d,:)/max(emp_phase_rate(d,:));
%     boot_phase_rate(d,:) = boot_phase_rate(d,:)/max(boot_phase_rate(d,:));
%     emp_tslu_rate(d,:) = emp_tslu_rate(d,:)/max(emp_tslu_rate(d,:));
   
boot_diff(d,:) = (emp_phase_rate(d,:)-boot_phase_rate(d,:))/mean_rate(d);

    
end
% errorbar(phase_hist,nanmean(emp_phase_rate),nanstd(emp_phase_rate)/sqrt(17))
% xlim([0 360])
% hold on
% errorbar(phase_hist,nanmean(boot_phase_rate),nanstd(boot_phase_rate)/sqrt(17),'r')
errorbar(phase_hist,nanmean(boot_diff),nanstd(boot_diff)/sqrt(17))