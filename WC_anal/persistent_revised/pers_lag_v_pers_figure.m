clear all
close all
cd C:\WC_Germany\persistent_revised

load corresponding_lfp_state_data

mec_cells = 1:21;
lec_cells = 22:28;
lec_cells(4) = []; %get rid of 12-09-A
% mec_cells(18) = []; %get rid of 04-07
mec_cells(21) = [];
for d = 1:28
    non_pers_states = find(mp_updur_corresp_lfp{d} < 1);
    mean_down_lag_np(d) = nanmean(lfp_down_lag{d}(non_pers_states));
end
[a,b] = corrcoef(mean_down_lag(mec_cells),pers_fract_across_cycles(mec_cells))
[a,b] = corrcoef(mean_down_lag(lec_cells),pers_fract_across_cycles(lec_cells))
[a,b] = corrcoef(mean_down_lag([mec_cells lec_cells]),pers_fract_across_cycles([mec_cells lec_cells]))
[a,b] = corrcoef(mean_down_lag_np(mec_cells),pers_fract_across_cycles(mec_cells))
[a,b] = corrcoef(mean_down_lag_np(lec_cells),pers_fract_across_cycles(lec_cells))
[a,b] = corrcoef(mean_down_lag_np([mec_cells lec_cells]),pers_fract_across_cycles([mec_cells lec_cells]))

plot(mean_down_lag(mec_cells),pers_fract_across_cycles(mec_cells),'o','markersize',8)
hold on
plot(mean_down_lag(lec_cells),pers_fract_across_cycles(lec_cells),'ro','markersize',8)
plot(mean_down_lag_np(mec_cells),pers_fract_across_cycles(mec_cells),'*','markersize',8)
plot(mean_down_lag_np(lec_cells),pers_fract_across_cycles(lec_cells),'r*','markersize',8)
xlim([0 1.5])
