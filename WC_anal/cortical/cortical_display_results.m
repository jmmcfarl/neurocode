clear all
close all

load C:\WC_Germany\Cortical_analysis\cortical_dir
load C:\WC_Germany\Cortical_analysis\sigmoid_fit\sigfitcompare
load C:\WC_Germany\Cortical_analysis\trig_avgs\trig_avg_data_sigfit
load C:\WC_Germany\Cortical_analysis\trig_avgs\trig_avg_data_wideband
load C:\WC_Germany\Cortical_analysis\trig_avgs\lf8_trig_avg_data_sigfit

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

l_5_pyr_fro = [l_5_pyr_fro l_5_pyr_pre];
l_23_pyr_fro = [l_23_pyr_fro l_23_pyr_pre];

norm_lf8_utrig_spk_t10 = lf8_utrig_spk_t10./repmat(max(lf8_utrig_spk_t10,[],2),1,length(lags));
norm_lf8_utrig_spk_t50 = lf8_utrig_spk_t50./repmat(max(lf8_utrig_spk_t50,[],2),1,length(lags));
norm_lf8_utrig_spk_t90 = lf8_utrig_spk_t90./repmat(max(lf8_utrig_spk_t90,[],2),1,length(lags));

%%
Fsd = 2016/8;
data_mat = norm_lf8_utrig_spk_t10;
xvec = lags/Fsd;

figure
[mean_vec8,se_vec8] = get_data_mean_se(data_mat,l_23_pyr_par);
h = errorbar(xvec,mean_vec8,se_vec8)
hold on
errorbar_tick(h,0.001,'units')

[mean_vec8,se_vec8] = get_data_mean_se(data_mat,l_56_pyr_par);
h= errorbar(xvec,mean_vec8,se_vec8,'r')
errorbar_tick(h,0.001,'units')

[mean_vec8,se_vec8] = get_data_mean_se(data_mat,l_23_pyr_fro);
h = errorbar(xvec,mean_vec8,se_vec8,'g')
errorbar_tick(h,0.001,'units')

[mean_vec8,se_vec8] = get_data_mean_se(data_mat,l_5_pyr_fro);
h = errorbar(xvec,mean_vec8,se_vec8,'k')
errorbar_tick(h,0.001,'units')

legend('L-23-par','L-56-par','L-23-Fro','L-5-Fro')
% xlim([-1 0.5])
grid
