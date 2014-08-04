clear all
close all

load C:\WC_Germany\Cortical_analysis\cortical_dir
load C:\WC_Germany\Cortical_analysis\trig_avgs\lf8_trig_avg_data_sigfit
load C:\WC_Germany\Cortical_analysis\cortical_mean_rate_data

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

view_all_cell_types(lf8_utrig_spk_t50,lags/Fsd)
line([0 0],[0 57],'Color','w')
xlim([-0.5 0.5])


norm_lf8_spk_t50 = lf8_utrig_spk_t50./repmat(max(lf8_utrig_spk_t50,[],2),1,1513);
norm_lf8_spk_t10 = lf8_utrig_spk_t10./repmat(max(lf8_utrig_spk_t10,[],2),1,1513);

cur_data = lf8_utrig_spk_t10;

figure
[mean_vec8,se_vec8] = get_data_mean_se(cur_data,l_23_pyr_par);
h = errorbar(lags/Fsd,mean_vec8,se_vec8)
hold on
errorbar_tick(h,0.001,'units')

[mean_vec8,se_vec8] = get_data_mean_se(cur_data,l_56_pyr_par);
h= errorbar(lags/Fsd,mean_vec8,se_vec8,'r')
errorbar_tick(h,0.001,'units')

[mean_vec8,se_vec8] = get_data_mean_se(cur_data,l_23_pyr_fro);
h = errorbar(lags/Fsd,mean_vec8,se_vec8,'g')
errorbar_tick(h,0.001,'units')

[mean_vec8,se_vec8] = get_data_mean_se(cur_data,l_5_pyr_fro);
h = errorbar(lags/Fsd,mean_vec8,se_vec8,'k')
errorbar_tick(h,0.001,'units')

[mean_vec8,se_vec8] = get_data_mean_se(cur_data,l_23_pyr_pre);
h = errorbar(lags/Fsd,mean_vec8,se_vec8,'c')
errorbar_tick(h,0.001,'units')

[mean_vec8,se_vec8] = get_data_mean_se(cur_data,l_5_pyr_pre);
h = errorbar(lags/Fsd,mean_vec8,se_vec8,'Color',[0.8 0.2 0])
errorbar_tick(h,0.001,'units')

legend('L-23-par','L-56-par','L-23-Fro','L-5-Fro','L-23-pre','L-5-pre')
xlim([-0.5 0.5])
grid

