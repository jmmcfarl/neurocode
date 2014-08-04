clear all
close all

load('C:\WC_Germany\overall_calcs\overall_dir.mat')
load C:\WC_Germany\overall_calcs\trig_avgs\trig_avg_data
load C:\WC_Germany\overall_calcs\trig_avgs\trig_avg_midpts

diff_types = unique(cell_type);

for i = 1:length(diff_types)
    
    mp_u_mp(i) = mean(lags(mp_utrig_mp_mid(cell_type==diff_types(i))));
    mp_u_mp_s(i) = std(lags(mp_utrig_mp_mid(cell_type==diff_types(i))));
    mp_u_lf8(i) = mean(lags(mp_utrig_lf8_mid(cell_type==diff_types(i))));
    mp_u_lf8_s(i) = std(lags(mp_utrig_lf8_mid(cell_type==diff_types(i))));
    mp_u_lf2(i) = mean(lags(mp_utrig_lf2_mid(cell_type==diff_types(i))));
    mp_u_lf2_s(i) = std(lags(mp_utrig_lf2_mid(cell_type==diff_types(i))));
    
    mp_d_mp(i) = mean(lags(mp_dtrig_mp_mid(cell_type==diff_types(i))));
    mp_d_mp_s(i) = std(lags(mp_dtrig_mp_mid(cell_type==diff_types(i))));
    mp_d_lf8(i) = mean(lags(mp_dtrig_lf8_mid(cell_type==diff_types(i))));
    mp_d_lf8_s(i) = std(lags(mp_dtrig_lf8_mid(cell_type==diff_types(i))));
    mp_d_lf2(i) = mean(lags(mp_dtrig_lf2_mid(cell_type==diff_types(i))));
    mp_d_lf2_s(i) = std(lags(mp_dtrig_lf2_mid(cell_type==diff_types(i))));

    lf8_u_mp(i) = mean(lags(lf8_utrig_mp_mid(cell_type==diff_types(i))));
    lf8_u_mp_s(i) = std(lags(lf8_utrig_mp_mid(cell_type==diff_types(i))));
    lf8_u_lf8(i) = mean(lags(lf8_utrig_lf8_mid(cell_type==diff_types(i))));
    lf8_u_lf8_s(i) = std(lags(lf8_utrig_lf8_mid(cell_type==diff_types(i))));
    lf8_u_lf2(i) = mean(lags(lf8_utrig_lf2_mid(cell_type==diff_types(i))));
    lf8_u_lf2_s(i) = std(lags(lf8_utrig_lf2_mid(cell_type==diff_types(i))));
    
    lf8_d_mp(i) = mean(lags(lf8_dtrig_mp_mid(cell_type==diff_types(i))));
    lf8_d_mp_s(i) = std(lags(lf8_dtrig_mp_mid(cell_type==diff_types(i))));
    lf8_d_lf8(i) = mean(lags(lf8_dtrig_lf8_mid(cell_type==diff_types(i))));
    lf8_d_lf8_s(i) = std(lags(lf8_dtrig_lf8_mid(cell_type==diff_types(i))));
    lf8_d_lf2(i) = mean(lags(lf8_dtrig_lf2_mid(cell_type==diff_types(i))));
    lf8_d_lf2_s(i) = std(lags(lf8_dtrig_lf2_mid(cell_type==diff_types(i))));

end

