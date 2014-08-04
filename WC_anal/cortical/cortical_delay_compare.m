clear all
close all
cd C:\WC_Germany\Cortical_analysis

load C:\WC_Germany\Cortical_analysis\trig_avgs\trig_avg_data_sigfit

Fsd = 2016/8;

l_23_pyr_par = 1:11;

int_par = 12:14;

l_56_pyr_par = 15:24;

l_23_pyr_fro = 25:28;

l_5_pyr_fro = 29:31;

l_23_pyr_pre = 32:36;

l_5_pyr_pre = 37:39;

l_23_pyr_fro = [l_23_pyr_fro l_23_pyr_pre];
l_5_pyr_fro = [l_5_pyr_fro l_5_pyr_pre];

thal = 40:43;

barrel = 44:47;

m_23_pyr_par_mp = nanmean(mp_utrig_mp_t10(l_23_pyr_par,:));
se_23_pyr_par_mp = nanstd(mp_utrig_mp_t10(l_23_pyr_par,:))/sqrt(length(l_23_pyr_par));

m_56_pyr_par_mp = nanmean(mp_utrig_mp_t10(l_56_pyr_par,:));
se_56_pyr_par_mp = nanstd(mp_utrig_mp_t10(l_56_pyr_par,:))/sqrt(length(l_56_pyr_par));

m_int_par_mp = nanmean(mp_utrig_mp_t10(int_par,:));
se_int_par_mp = nanstd(mp_utrig_mp_t10(int_par,:))/sqrt(length(int_par));

m_l_23_pyr_fro_mp = nanmean(mp_utrig_mp_t10(l_23_pyr_fro,:));
se_l_23_pyr_fro_mp = nanstd(mp_utrig_mp_t10(l_23_pyr_fro,:))/sqrt(length(l_23_pyr_fro));

m_l_5_pyr_fro_mp = nanmean(mp_utrig_mp_t10(l_5_pyr_fro,:));
se_l_5_pyr_fro_mp = nanstd(mp_utrig_mp_t10(l_5_pyr_fro,:))/sqrt(length(l_5_pyr_fro));

m_l_23_pyr_pre_mp = nanmean(mp_utrig_mp_t10(l_23_pyr_pre,:));
se_l_23_pyr_pre_mp = nanstd(mp_utrig_mp_t10(l_23_pyr_pre,:))/sqrt(length(l_23_pyr_pre));

m_l_5_pyr_pre_mp = nanmean(mp_utrig_mp_t10(l_5_pyr_pre,:));
se_l_5_pyr_pre_mp = nanstd(mp_utrig_mp_t10(l_5_pyr_pre,:))/sqrt(length(l_5_pyr_pre));

m_thal_mp = nanmean(mp_utrig_mp_t10(thal,:));
se_thal_mp = nanstd(mp_utrig_mp_t10(thal,:))/sqrt(length(thal));

m_barrel_mp = nanmean(mp_utrig_mp_t10(barrel,:));
se_barrel_mp = nanstd(mp_utrig_mp_t10(barrel,:))/sqrt(length(barrel));



m_23_pyr_par_lf8 = nanmean(mp_utrig_lf8_t10(l_23_pyr_par,:));
se_23_pyr_par_lf8 = nanstd(mp_utrig_lf8_t10(l_23_pyr_par,:))/sqrt(length(l_23_pyr_par));

m_56_pyr_par_lf8 = nanmean(mp_utrig_lf8_t10(l_56_pyr_par,:));
se_56_pyr_par_lf8 = nanstd(mp_utrig_lf8_t10(l_56_pyr_par,:))/sqrt(length(l_56_pyr_par));

m_int_par_lf8 = nanmean(mp_utrig_lf8_t10(int_par,:));
se_int_par_lf8 = nanstd(mp_utrig_lf8_t10(int_par,:))/sqrt(length(int_par));

m_l_23_pyr_fro_lf8 = nanmean(mp_utrig_lf8_t10(l_23_pyr_fro,:));
se_l_23_pyr_fro_lf8 = nanstd(mp_utrig_lf8_t10(l_23_pyr_fro,:))/sqrt(length(l_23_pyr_fro));

m_l_5_pyr_fro_lf8 = nanmean(mp_utrig_lf8_t10(l_5_pyr_fro,:));
se_l_5_pyr_fro_lf8 = nanstd(mp_utrig_lf8_t10(l_5_pyr_fro,:))/sqrt(length(l_5_pyr_fro));

m_l_23_pyr_pre_lf8 = nanmean(mp_utrig_lf8_t10(l_23_pyr_pre,:));
se_l_23_pyr_pre_lf8 = nanstd(mp_utrig_lf8_t10(l_23_pyr_pre,:))/sqrt(length(l_23_pyr_pre));

m_l_5_pyr_pre_lf8 = nanmean(mp_utrig_lf8_t10(l_5_pyr_pre,:));
se_l_5_pyr_pre_lf8 = nanstd(mp_utrig_lf8_t10(l_5_pyr_pre,:))/sqrt(length(l_5_pyr_pre));

m_thal_lf8 = nanmean(mp_utrig_lf8_t10(thal,:));
se_thal_lf8 = nanstd(mp_utrig_lf8_t10(thal,:))/sqrt(length(thal));

m_barrel_lf8 = nanmean(mp_utrig_lf8_t10(barrel,:));
se_barrel_lf8 = nanstd(mp_utrig_lf8_t10(barrel,:))/sqrt(length(barrel));


m_23_pyr_par_lf5 = nanmean(mp_utrig_lf5_t10(l_23_pyr_par,:));
se_23_pyr_par_lf5 = nanstd(mp_utrig_lf5_t10(l_23_pyr_par,:))/sqrt(length(l_23_pyr_par));

m_56_pyr_par_lf5 = nanmean(mp_utrig_lf5_t10(l_56_pyr_par,:));
se_56_pyr_par_lf5 = nanstd(mp_utrig_lf5_t10(l_56_pyr_par,:))/sqrt(length(l_56_pyr_par));

m_int_par_lf5 = nanmean(mp_utrig_lf5_t10(int_par,:));
se_int_par_lf5 = nanstd(mp_utrig_lf5_t10(int_par,:))/sqrt(length(int_par));

m_l_23_pyr_fro_lf5 = nanmean(mp_utrig_lf5_t10(l_23_pyr_fro,:));
se_l_23_pyr_fro_lf5 = nanstd(mp_utrig_lf5_t10(l_23_pyr_fro,:))/sqrt(length(l_23_pyr_fro));

m_l_5_pyr_fro_lf5 = nanmean(mp_utrig_lf5_t10(l_5_pyr_fro,:));
se_l_5_pyr_fro_lf5 = nanstd(mp_utrig_lf5_t10(l_5_pyr_fro,:))/sqrt(length(l_5_pyr_fro));

m_l_23_pyr_pre_lf5 = nanmean(mp_utrig_lf5_t10(l_23_pyr_pre,:));
se_l_23_pyr_pre_lf5 = nanstd(mp_utrig_lf5_t10(l_23_pyr_pre,:))/sqrt(length(l_23_pyr_pre));

m_l_5_pyr_pre_lf5 = nanmean(mp_utrig_lf5_t10(l_5_pyr_pre,:));
se_l_5_pyr_pre_lf5 = nanstd(mp_utrig_lf5_t10(l_5_pyr_pre,:))/sqrt(length(l_5_pyr_pre));

m_thal_lf5 = nanmean(mp_utrig_lf5_t10(thal,:));
se_thal_lf5 = nanstd(mp_utrig_lf5_t10(thal,:))/sqrt(length(thal));

m_barrel_lf5 = nanmean(mp_utrig_lf5_t10(barrel,:));
se_barrel_lf5 = nanstd(mp_utrig_lf5_t10(barrel,:))/sqrt(length(barrel));



m_23_pyr_par_spk = nanmean(mp_utrig_spk_t10(l_23_pyr_par,:));
se_23_pyr_par_spk = nanstd(mp_utrig_spk_t10(l_23_pyr_par,:))/sqrt(length(l_23_pyr_par));

m_56_pyr_par_spk = nanmean(mp_utrig_spk_t10(l_56_pyr_par,:));
se_56_pyr_par_spk = nanstd(mp_utrig_spk_t10(l_56_pyr_par,:))/sqrt(length(l_56_pyr_par));

m_int_par_spk = nanmean(mp_utrig_spk_t10(int_par,:));
se_int_par_spk = nanstd(mp_utrig_spk_t10(int_par,:))/sqrt(length(int_par));

m_l_23_pyr_fro_spk = nanmean(mp_utrig_spk_t10(l_23_pyr_fro,:));
se_l_23_pyr_fro_spk = nanstd(mp_utrig_spk_t10(l_23_pyr_fro,:))/sqrt(length(l_23_pyr_fro));

m_l_5_pyr_fro_spk = nanmean(mp_utrig_spk_t10(l_5_pyr_fro,:));
se_l_5_pyr_fro_spk = nanstd(mp_utrig_spk_t10(l_5_pyr_fro,:))/sqrt(length(l_5_pyr_fro));

m_l_23_pyr_pre_spk = nanmean(mp_utrig_spk_t10(l_23_pyr_pre,:));
se_l_23_pyr_pre_spk = nanstd(mp_utrig_spk_t10(l_23_pyr_pre,:))/sqrt(length(l_23_pyr_pre));

m_l_5_pyr_pre_spk = nanmean(mp_utrig_spk_t10(l_5_pyr_pre,:));
se_l_5_pyr_pre_spk = nanstd(mp_utrig_spk_t10(l_5_pyr_pre,:))/sqrt(length(l_5_pyr_pre));

m_thal_spk = nanmean(mp_utrig_spk_t10(thal,:));
se_thal_spk = nanstd(mp_utrig_spk_t10(thal,:))/sqrt(length(thal));

m_barrel_spk = nanmean(mp_utrig_spk_t10(barrel,:));
se_barrel_spk = nanstd(mp_utrig_spk_t10(barrel,:))/sqrt(length(barrel));


figure
title('T10 triggered MP')
plot(lags/Fsd,m_23_pyr_par_mp)
hold on
plot(lags/Fsd,m_56_pyr_par_mp,'k')
plot(lags/Fsd,m_int_par_mp,'r')
plot(lags/Fsd,m_l_23_pyr_fro_mp,'c')
plot(lags/Fsd,m_l_5_pyr_fro_mp,'g')
plot(lags/Fsd,m_thal_mp,'Color',[0.3 0.3 0.5])
plot(lags/Fsd,m_barrel_mp,'y')
legend('23 Parietal','56 Parietal','Parietal Int','23 frontal','5 frontal','thalamic','barrel')
xlim([-0.5 0.5])
grid

figure
title('T10 triggered LF8')
plot(lags/Fsd,m_23_pyr_par_lf8)
hold on
plot(lags/Fsd,m_56_pyr_par_lf8,'k')
plot(lags/Fsd,m_int_par_lf8,'r')
plot(lags/Fsd,m_l_23_pyr_fro_lf8,'c')
plot(lags/Fsd,m_l_5_pyr_fro_lf8,'g')
plot(lags/Fsd,m_thal_lf8,'Color',[0.3 0.3 0.5])
plot(lags/Fsd,m_barrel_lf8,'y')
legend('23 Parietal','56 Parietal','Parietal Int','23 frontal','5 frontal','thalamic','barrel')
xlim([-0.5 0.5])
grid

figure
title('T10 triggered Spk')
plot(lags/Fsd,m_23_pyr_par_spk)
hold on
plot(lags/Fsd,m_56_pyr_par_spk,'k')
plot(lags/Fsd,m_int_par_spk,'r')
plot(lags/Fsd,m_l_23_pyr_fro_spk,'c')
plot(lags/Fsd,m_l_5_pyr_fro_spk,'g')
plot(lags/Fsd,m_thal_spk,'Color',[0.3 0.3 0.5])
plot(lags/Fsd,m_barrel_spk,'y')
legend('23 Parietal','56 Parietal','Parietal Int','23 frontal','5 frontal','thalamic','barrel')
xlim([-0.5 0.5])
grid



