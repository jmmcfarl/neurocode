clear all
close all
cd C:\WC_Germany\Cortical_analysis

load C:\WC_Germany\Cortical_analysis\cortical_dir
load C:\WC_Germany\Cortical_analysis\UDS_dur_raw\UDS_raw_data
load C:\WC_Germany\Cortical_analysis\UDS_synch_state_dur\UDS_synch_state_dur_data
load C:\WC_Germany\Cortical_analysis\lf8_period_f

dsf = 8;
Fsd = 2016/dsf;


for d = 1:length(sess_data)
    d
    lfp_up_lag_near{d} = zeros(size(synch_ups{d}));
    lfp_up_lag_prev{d} = zeros(size(synch_ups{d}));
    lfp_up_lag_near5{d} = zeros(size(synch_ups{d}));
    lfp_up_lag_prev5{d} = zeros(size(synch_ups{d}));

    lfp_down_lag_near{d} = zeros(size(synch_ups{d}));
    lfp_down_lag_prev{d} = zeros(size(synch_ups{d}));
    mp_up_phase{d} = zeros(size(synch_ups{d}));
    mp_down_phase{d} = zeros(size(synch_ups{d}));
    mp_updur_lfp_cor{d} = zeros(size(synch_ups{d}));
    mp_updur_lfp_uncor{d} = zeros(size(synch_ups{d}));

    for i = 1:length(synch_ups{d})
    
        mp_up_phase{d}(i) = lf8_period_p{d}(up_trans{d}(synch_ups{d}(i)));
        mp_down_phase{d}(i) = lf8_period_p{d}(down_trans{d}(synch_ups{d}(i)));
        
        %acquire MP up related lfp transitions
        [dummy,near_lfp_up] = min(abs(up_trans{d}(synch_ups{d}(i)) - up_trans8{d}));
        [dummy,near_lfp_down] = min(abs(down_trans{d}(synch_ups{d}(i)) - down_trans8{d}));
        [dummy,near_lfp_up5] = min(abs(up_trans{d}(synch_ups{d}(i)) - up_trans5{d}));
        [dummy,near_lfp_down5] = min(abs(down_trans{d}(synch_ups{d}(i)) - down_trans5{d}));

        if ismember(near_lfp_up,synch_ups8{d})
            lfp_up_lag_near{d}(i) = (up_trans{d}(synch_ups{d}(i))-up_trans8{d}(near_lfp_up))/Fsd;
        else
            lfp_up_lag_near{d}(i) = nan;
        end
         if ismember(near_lfp_up5,synch_ups8{d})
            lfp_up_lag_near5{d}(i) = (up_trans{d}(synch_ups{d}(i))-up_trans5{d}(near_lfp_up5))/Fsd;
        else
            lfp_up_lag_near5{d}(i) = nan;
        end
       
        if ismember(near_lfp_down,synch_downs8{d})
            lfp_down_lag_near{d}(i) = (down_trans{d}(synch_ups{d}(i))-down_trans8{d}(near_lfp_down))/Fsd; 
        else
            lfp_down_lag_near{d}(i) = nan;
        end
         if ismember(near_lfp_down5,synch_downs8{d})
            lfp_down_lag_near5{d}(i) = (down_trans{d}(synch_ups{d}(i))-down_trans5{d}(near_lfp_down5))/Fsd; 
        else
            lfp_down_lag_near5{d}(i) = nan;
        end
       
        prev_lfp_up = find(up_trans8{d}<up_trans{d}(synch_ups{d}(i)),1,'last');
       
        if ~isempty(prev_lfp_up) & ismember(prev_lfp_up,synch_ups8{d})
            
            lfp_up_lag_prev{d}(i) = (up_trans{d}(synch_ups{d}(i))-up_trans8{d}(prev_lfp_up))/Fsd;
            mp_updur_lfp_cor{d}(i) = lf8_period_f{d}(down_trans{d}(synch_ups{d}(i))-lfp_up_lag_near{d}(i)*Fsd)...
                -lf8_period_f{d}(up_trans{d}(synch_ups{d}(i))-lfp_up_lag_near{d}(i)*Fsd);
            mp_updur_lfp_uncor{d}(i) = lf8_period_f{d}(down_trans{d}(synch_ups{d}(i)))...
                -lf8_period_f{d}(up_trans{d}(synch_ups{d}(i)));
            
        else
            
            lfp_up_lag_prev{d}(i) = nan;
            mp_updur_lfp_cor{d}(i) = nan;
            mp_updur_lfp_uncor{d}(i) = nan;
            
        end

        
        


    end


    pers_fract_cor(d) = length(find(mp_updur_lfp_cor{d}>1))/length(mp_updur_lfp_cor{d});
    pers_fract_uncor(d) = length(find(mp_updur_lfp_uncor{d}>1))/length(mp_updur_lfp_uncor{d});


end

save C:\WC_Germany\Cortical_analysis\lag_phase_data pers_fract* mp_updur* lfp_up_lag* lfp_down_lag* mp*phase
