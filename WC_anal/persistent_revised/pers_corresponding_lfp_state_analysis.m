clear all

load C:\WC_Germany\persistent_revised\pers_revised_dir
load C:\WC_Germany\persistent_revised\UDS_dur_raw\UDS_raw_data
load C:\WC_Germany\persistent_revised\UDS_synch_state_dur\UDS_synch_state_dur_data
load C:\WC_Germany\persistent_revised\lf8_period_f
Fsd = 2016/8;

for d = 1:28
    d
    lfp_up_lag{d} = zeros(size(synch_ups{d}));
    lfp_down_lag{d} = zeros(size(synch_ups{d}));
    mp_updur_corresp_lfp{d} = zeros(size(synch_ups{d}));

    for i = 1:length(synch_ups{d})

        %acquire MP up related lfp transitions
        [dummy,near_lfp_up] = min(abs(up_trans{d}(synch_ups{d}(i)) - up_trans8{d}));

        if ismember(near_lfp_up,synch_ups8{d})
            lfp_up_lag{d}(i) = (up_trans{d}(synch_ups{d}(i))-up_trans8{d}(near_lfp_up))/Fsd;
            mp_updur_corresp_lfp{d}(i) = lf8_period_f{d}(down_trans{d}(synch_ups{d}(i))-lfp_up_lag{d}(i)*Fsd)...
                - lf8_period_f{d}(up_trans{d}(synch_ups{d}(i))-lfp_up_lag{d}(i)*Fsd);
            %find preceding lfp down transition after the corresponding up
            %transition
            corresp_lfp_down = find(down_trans8{d} > up_trans8{d}(near_lfp_up) ...
                & down_trans8{d} < down_trans{d}(synch_ups{d}(i)),1,'last');
            if ~isempty(corresp_lfp_down)
                lfp_down_lag{d}(i) = (down_trans{d}(synch_ups{d}(i)) - down_trans8{d}(corresp_lfp_down))/Fsd;
            else
                lfp_down_lag{d}(i) = nan;
            end
        else
            lfp_up_lag{d}(i) = nan;
            lfp_down_lag{d}(i) = nan;
            mp_updur_corresp_lfp{d}(i) = nan;
        end

    end

    pers_fract_within_cycle(d) = length(find(mp_updur_corresp_lfp{d}>0.5))/length(mp_updur_corresp_lfp{d});
    pers_fract_across_cycles(d) = length(find(mp_updur_corresp_lfp{d} > 1))/length(mp_updur_corresp_lfp{d});
    pers_fract_within_np(d) = length(find(mp_updur_corresp_lfp{d} > 0.5 & mp_updur_corresp_lfp{d} < 1))...
        /length(find(mp_updur_corresp_lfp{d} < 1));
    persistent_dur(d) = nanmean(up_state_dur{d}(synch_ups{d}(mp_updur_corresp_lfp{d} > 1)));
    mean_up_lag(d) = nanmean(lfp_up_lag{d});
    mean_down_lag(d) = nanmean(lfp_down_lag{d});
    median_up_lag(d) = nanmedian(lfp_up_lag{d});
    median_down_lag(d) = nanmedian(lfp_down_lag{d});
end

mec_cells = 1:20;
lec_cells = 22:28;
lec_cells(4) = []; %get rid of 12-09-A

save corresponding_lfp_state_data pers_fract* mean*lag median*lag lfp*lag mp_updur*