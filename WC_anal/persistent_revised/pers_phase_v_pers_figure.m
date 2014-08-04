clear all
close all
cd C:\WC_Germany\persistent_revised

load lag_phase_data
for i = 1:length(mp_up_phase)
    std_up_phase(i) = 180/pi*circ_std(mp_up_phase{i}*2*pi/360);
    std_down_phase(i) = 180/pi*circ_std(mp_down_phase{i}*2*pi/360);
    mean_up_phase(i) = 180/pi*circ_mean(mp_up_phase{i}*2*pi/360);
    mean_down_phase(i) = 180/pi*circ_mean(mp_down_phase{i}*2*pi/360);
end

mec_cells = 1:21;
lec_cells = 22:28;
lec_cells(4) = []; %get rid of 12-09-A
% mec_cells(18) = []; %get rid of 04-07
mec_cells(21) = [];

% mean_up_phase = mean_up_phase([mec_cells lec_cells]);
% mean_down_phase = mean_down_phase([mec_cells lec_cells]);
% pers_fract_uncor = pers_fract_uncor([mec_cells lec_cells]);

[rho_mec pval_mec] = circ_corrcl(mean_down_phase(mec_cells)/360*2*pi, pers_fract_uncor(mec_cells))
[rho_lec pval_lec] = circ_corrcl(mean_down_phase(lec_cells)/360*2*pi, pers_fract_uncor(lec_cells))
[rho_all pval_all] = circ_corrcl(mean_down_phase([mec_cells lec_cells])/360*2*pi, pers_fract_uncor([mec_cells lec_cells]))

mean_down_phase(mean_down_phase < 0) = mean_down_phase(mean_down_phase < 0) + 360;
plot(mean_down_phase(mec_cells),pers_fract_uncor(mec_cells),'o','markersize',8)
hold on
plot(mean_down_phase(lec_cells),pers_fract_uncor(lec_cells),'ro','markersize',8)


m_mec_up_phase = 180/pi*circ_mean(mean_up_phase(mec_cells)*pi/180);
s_mec_up_phase = 180/pi*circ_std(mean_up_phase(mec_cells)*pi/180)/sqrt(length(mec_cells));

m_lec_up_phase = 180/pi*circ_mean(mean_up_phase(lec_cells)*pi/180);
s_lec_up_phase = 180/pi*circ_std(mean_up_phase(lec_cells)*pi/180)/sqrt(length(lec_cells));

m_mec_down_phase = 180/pi*circ_mean(mean_down_phase(mec_cells)*pi/180);
s_mec_down_phase = 180/pi*circ_std(mean_down_phase(mec_cells)*pi/180)/sqrt(length(mec_cells));

m_lec_down_phase = 180/pi*circ_mean(mean_down_phase(lec_cells)*pi/180);
s_lec_down_phase = 180/pi*circ_std(mean_down_phase(lec_cells)*pi/180)/sqrt(length(lec_cells));