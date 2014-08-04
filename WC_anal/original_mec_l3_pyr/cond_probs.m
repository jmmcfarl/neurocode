%% calculate conditional up and down probabilities
clear all
% load C:\WC_Germany\JMM_Analysis_pyr\UDS_dur_run_hist_v2\up_per_data_10_28

load C:\WC_Germany\JMM_Analysis_pyr\UDS_dur_run_hist_v2\UDS_dur_data_over_smooth
load C:\WC_Germany\JMM_Analysis_pyr\dir_tree_update

dsf = 8;
Fsd = 2016/dsf;
lfp_lags = 0; 
lfp_lags = round(lfp_lags*Fsd);

for d = 1:length(dir_array)
    
    cd(dir_array{d})
    pwd
    load used_data wcv_minus_spike
    
    wcv_d = downsample(wcv_minus_spike,dsf);
    
    for l = 1:length(lfp_lags)
    
    wcv_up = zeros(length(wcv_d),1);
    lf8_up = wcv_up;
    temp_lf8_up_trans = up_trans8{d};
    temp_lf8_down_trans = down_trans8{d};
    temp_lf8_up_trans = temp_lf8_up_trans - lfp_lags(l);
    temp_lf8_down_trans = temp_lf8_down_trans-lfp_lags(l);
    
      %find wcv up times
    for i = 1:length(up_trans{d})
    
        wcv_up(up_trans{d}(i):down_trans{d}(i)) = 1;
        
    end
    
    %find lf8 up times
    for i = 1:length(temp_lf8_up_trans)
       
        lf8_up(temp_lf8_up_trans(i):temp_lf8_down_trans(i)) =1;
        
    end
    
    wcv_up_prob(d) = sum(wcv_up)/length(wcv_up);
    wcv_down_prob(d) = 1-wcv_up_prob(d);
    lf8_up_prob(d) = sum(lf8_up)/length(lf8_up);
    lf8_down_prob(d) = 1-lf8_up_prob(d);
    
    wup_ldown = wcv_up & ~lf8_up;
    wup_lup = wcv_up & lf8_up;
    wdown_lup = ~wcv_up & lf8_up;
    wdown_ldown = ~wcv_up & ~lf8_up;
    
    wup_cond_lup(l,d) = sum(wup_lup)/sum(lf8_up);
    wdown_cond_lup(l,d) = 1-wup_cond_lup(l,d);
    wup_cond_ldown(l,d) = sum(wup_ldown)/sum(~lf8_up);
    wdown_cond_ldown(l,d) = 1-wup_cond_ldown(l,d);
    lup_cond_wup(l,d) = sum(wup_lup)/sum(wcv_up);
    ldown_cond_wup(l,d) = 1-lup_cond_wup(l,d);
    lup_cond_wdown(l,d) = sum(wdown_lup)/sum(~wcv_up);
    ldown_cond_wdown(l,d) = 1-lup_cond_wdown(l,d);
   
    end
    
end

cd C:\WC_Germany\JMM_Analysis_pyr