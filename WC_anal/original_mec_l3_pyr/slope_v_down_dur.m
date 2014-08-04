clear all

load C:\WC_Germany\JMM_analysis_pyr\dir_tree_update
load C:\WC_Germany\JMM_analysis_pyr\UDS_dur_run_hist_v2\up_per_data_10_28
load C:\WC_Germany\JMM_analysis_pyr\t_90_data
load C:\WC_Germany\JMM_analysis_pyr\wcv_up_tau


for d = 1:length(dir_array)
    
    cur_down_durs = down_state_dur{d};
    cur_slopes = rltau_wup{d};
    cur_up_trans = up_trans{d};
    
   if ~exist('spike_time.mat')
        load spike_time_br
   else
       load spike_time
   end
    spike_id = round(spkid/8);
    isis = [0;diff(spike_id)];

    spike_lag = zeros(length(cur_up_trans),1);
    
    for i = 1:length(cur_up_trans)
       
        next_spike = find(spike_id > cur_up_trans(i),1,'first');
        if ~isempty(next_spike)
        spike_lag(i) = isis(next_spike);
        else
            spike_lag(i) = nan;
        end
        
    end
    
    good_pts = ~isnan(spike_lag') & ~isnan(cur_slopes);
    
    [a,b] = corrcoef(cur_slopes(good_pts),log(spike_lag(good_pts)));
    
    plot(log(spike_lag),cur_slopes,'o')
    title(sprintf('Corr = %0.2g   P = %0.2g',a(2,1),b(2,1)))
    t_names = ['C:\WC_Germany\JMM_analysis_pyr\State_Trig\isi_v_slope_' f_names{d}];
    print('-dpng',t_names)
    close all
    
    
    
    
    
    
    
end