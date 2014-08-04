clear all
close all

load C:\WC_Germany\Cortical_analysis\cortical_dir
% load C:\WC_Germany\Cortical_analysis\trig_avgs\trig_avg_data_wideband
load C:\WC_Germany\Cortical_analysis\trig_avgs\lf8_utrig_avg_spk
[dummy,layer_23] = find_sess_data_fields(sess_data,'layer','23');
[dummy,layer_56] = find_sess_data_fields(sess_data,'layer','56');
[dummy,layer_5] = find_sess_data_fields(sess_data,'layer','5');
[dummy,parietal] = find_sess_data_fields(sess_data,'region','parietal');
[dummy,prefrontal] = find_sess_data_fields(sess_data,'region','prefrontal');
[dummy,frontal] = find_sess_data_fields(sess_data,'region','frontal');

l_23_pyr_par = find(layer_23 & parietal);
l_56_pyr_par = find(layer_56 & parietal);
l_23_pyr_fro = find(layer_23 & frontal);
l_5_pyr_fro = find(layer_5 & frontal);
l_23_pyr_pre = find(layer_23 & prefrontal);
l_5_pyr_pre = find(layer_5 & prefrontal);


Fsd = 2016/8;

figure
[mean_vec8,se_vec8] = get_data_mean_se(mp_utrig_lf8,l_23_pyr_par);
h = errorbar(lags/Fsd,mean_vec8,se_vec8)
hold on
errorbar_tick(h,0.001,'units')

[mean_vec8,se_vec8] = get_data_mean_se(mp_utrig_lf8,l_56_pyr_par);
h= errorbar(lags/Fsd,mean_vec8,se_vec8,'r')
errorbar_tick(h,0.001,'units')

[mean_vec8,se_vec8] = get_data_mean_se(mp_utrig_lf8,l_23_pyr_fro);
h = errorbar(lags/Fsd,mean_vec8,se_vec8,'g')
errorbar_tick(h,0.001,'units')

[mean_vec8,se_vec8] = get_data_mean_se(mp_utrig_lf8,l_5_pyr_fro);
h = errorbar(lags/Fsd,mean_vec8,se_vec8,'k')
errorbar_tick(h,0.001,'units')

[mean_vec8,se_vec8] = get_data_mean_se(mp_utrig_lf8,l_23_pyr_pre);
h = errorbar(lags/Fsd,mean_vec8,se_vec8,'c')
errorbar_tick(h,0.001,'units')

[mean_vec8,se_vec8] = get_data_mean_se(mp_utrig_lf8,l_5_pyr_pre);
h = errorbar(lags/Fsd,mean_vec8,se_vec8,'Color',[0.8 0.2 0])
errorbar_tick(h,0.001,'units')

legend('L-23-par','L-56-par','L-23-Fro','L-5-Fro','L-23-pre','L-5-pre')
xlim([-1 1])
grid


figure
[mean_vec8,se_vec8] = get_data_mean_se(mp_utrig_lf5,l_23_pyr_par);
h = errorbar(lags/Fsd,mean_vec8,se_vec8)
hold on
errorbar_tick(h,0.001,'units')

[mean_vec8,se_vec8] = get_data_mean_se(mp_utrig_lf5,l_56_pyr_par);
h= errorbar(lags/Fsd,mean_vec8,se_vec8,'r')
errorbar_tick(h,0.001,'units')

[mean_vec8,se_vec8] = get_data_mean_se(mp_utrig_lf5,l_23_pyr_fro);
h = errorbar(lags/Fsd,mean_vec8,se_vec8,'g')
errorbar_tick(h,0.001,'units')

[mean_vec8,se_vec8] = get_data_mean_se(mp_utrig_lf5,l_5_pyr_fro);
h = errorbar(lags/Fsd,mean_vec8,se_vec8,'k')
errorbar_tick(h,0.001,'units')

[mean_vec8,se_vec8] = get_data_mean_se(mp_utrig_lf5,l_23_pyr_pre);
h = errorbar(lags/Fsd,mean_vec8,se_vec8,'c')
errorbar_tick(h,0.001,'units')

[mean_vec8,se_vec8] = get_data_mean_se(mp_utrig_lf5,l_5_pyr_pre);
h = errorbar(lags/Fsd,mean_vec8,se_vec8,'Color',[0.8 0.2 0])
errorbar_tick(h,0.001,'units')

legend('L-23-par','L-56-par','L-23-Fro','L-5-Fro','L-23-pre','L-5-pre')
xlim([-1 1])
grid


figure
[mean_vec8,se_vec8] = get_data_mean_se(mp_dtrig_lf8,l_23_pyr_par);
h = errorbar(lags/Fsd,mean_vec8,se_vec8)
hold on
errorbar_tick(h,0.001,'units')

[mean_vec8,se_vec8] = get_data_mean_se(mp_dtrig_lf8,l_56_pyr_par);
h= errorbar(lags/Fsd,mean_vec8,se_vec8,'r')
errorbar_tick(h,0.001,'units')

[mean_vec8,se_vec8] = get_data_mean_se(mp_dtrig_lf8,l_23_pyr_fro);
h = errorbar(lags/Fsd,mean_vec8,se_vec8,'g')
errorbar_tick(h,0.001,'units')

[mean_vec8,se_vec8] = get_data_mean_se(mp_dtrig_lf8,l_5_pyr_fro);
h = errorbar(lags/Fsd,mean_vec8,se_vec8,'k')
errorbar_tick(h,0.001,'units')

[mean_vec8,se_vec8] = get_data_mean_se(mp_dtrig_lf8,l_23_pyr_pre);
h = errorbar(lags/Fsd,mean_vec8,se_vec8,'c')
errorbar_tick(h,0.001,'units')

[mean_vec8,se_vec8] = get_data_mean_se(mp_dtrig_lf8,l_5_pyr_pre);
h = errorbar(lags/Fsd,mean_vec8,se_vec8,'Color',[0.8 0.2 0])
errorbar_tick(h,0.001,'units')

legend('L-23-par','L-56-par','L-23-Fro','L-5-Fro','L-23-pre','L-5-pre')
xlim([-1 1])
grid
