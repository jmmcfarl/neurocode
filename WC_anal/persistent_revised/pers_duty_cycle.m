clear all

load C:\WC_Germany\persistent_revised\pers_revised_dir
load C:\WC_Germany\persistent_revised\UDS_dur_raw\UDS_raw_data
load C:\WC_Germany\persistent_revised\UDS_synch_state_dur\UDS_synch_state_dur_data

for d = 1:28
    lf8_dc = zeros(length(synch_ups8{d})-1,1);
    for i = 1:length(synch_ups8{d})-1
        
        cur_up_dur = down_trans8{d}(synch_ups8{d}(i))- up_trans8{d}(synch_ups8{d}(i));
        cur_down_dur = up_trans8{d}(synch_ups8{d}(i+1))-down_trans8{d}(synch_ups8{d}(i+1)-1);
        lf8_dc(i) = cur_up_dur/(cur_up_dur+cur_down_dur);
    end
    
            synch_ups8{d}(synch_ups8{d} == length(up_state_dur8{d})) = [];
        overall_dc = up_state_dur8{d}(1:end-1)./(up_state_dur8{d}(1:end-1)+down_state_dur8{d});
        desynch_ups = setdiff(1:length(up_state_dur8{d}(1:end-1)),synch_ups8{d});
        desynch_downs = setdiff(1:length(down_state_dur8{d}),synch_downs8{d});
        either_desynch = unique([desynch_ups desynch_downs]);
        overall_dc(either_desynch) = [];

    
    lf8_duty_cycle(d)=  nanmean(lf8_dc);
    lf8_duty_cycle2(d) = nanmean(overall_dc);
    
end

mec_cells = 1:21;
lec_cells = 22:28;
% lec_cells(4) = []; %get rid of 12-09-A
% mec_cells(18) = []; %get rid of 04-07

lf8_duty_cycle = lf8_duty_cycle([mec_cells lec_cells]);
lf8_duty_cycle2 = lf8_duty_cycle2([mec_cells lec_cells]);