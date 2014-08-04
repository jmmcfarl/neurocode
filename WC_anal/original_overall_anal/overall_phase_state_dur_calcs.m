clear all

load('C:\WC_Germany\overall_calcs\overall_dir.mat')
load C:\WC_Germany\overall_calcs\UDS_dur_raw\UDS_raw_data
load C:\WC_Germany\overall_calcs\lf8_period_f_data
load C:\WC_Germany\overall_calcs\UDS_synch_state_dur\UDS_synch_state_dur_data

dsf = 8;
Fsd = 2016/dsf;


for d = 1:length(over_dir)
    d
    lfp_up_lag_near{d} = zeros(size(synch_ups{d}));
    lfp_up_lag_prev{d} = zeros(size(synch_ups{d}));
    mp_updur_lfp_cor{d} = zeros(size(synch_ups{d}));
        mp_updur_lfp_uncor{d} = zeros(size(synch_ups{d}));
    mp_up_phase{d} = zeros(size(synch_ups{d}));
    mp_down_phase{d} = zeros(size(synch_ups{d}));
    
    for i = 1:length(synch_ups{d})

        
        mp_up_phase{d}(i) = lf8_period_p{d}(up_trans{d}(synch_ups{d}(i)));
        mp_down_phase{d}(i) = lf8_period_p{d}(down_trans{d}(synch_ups{d}(i)));
        
        %acquire MP up related lfp transitions
        [dummy,near_lfp_up] = min(abs(up_trans{d}(synch_ups{d}(i)) - up_trans8{d}(synch_ups8{d})));
                [dummy,near_lfp_down] = min(abs(down_trans{d}(synch_ups{d}(i)) - down_trans8{d}(synch_ups8{d})));

        lfp_up_lag_near{d}(i) = (up_trans{d}(synch_ups{d}(i))-up_trans8{d}(synch_ups8{d}(near_lfp_up)))/Fsd;
        lfp_down_lag_near{d}(i) = (down_trans{d}(synch_ups{d}(i))-down_trans8{d}(synch_ups8{d}(near_lfp_down)))/Fsd;

        prev_lfp_up = find(up_trans8{d}<up_trans{d}(synch_ups{d}(i)),1,'last');
       
        if ~isempty(prev_lfp_up)
            
            lfp_up_lag_prev{d}(i) = (up_trans{d}(i)-up_trans8{d}(prev_lfp_up))/Fsd;
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

save C:\WC_Germany\overall_calcs\lag_phase_data mp_updur_lfp* pers_fract* lfp_up_lag* lfp_down_lag* mp*phase
