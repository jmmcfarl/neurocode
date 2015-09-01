clear all
close all
%%
load F:\WC_Germany\overall_EC\overall_allcells_dir
addpath('F:\Code\WC_anal\general\')
addpath('F:\WC_Germany\Overall_EC\')
addpath('F:\Code\Chronux\spectral_analysis\continuous\')
addpath('F:\WC_Germany\hsmm_state_detection\\')
addpath('F:\WC_Germany\parietal_cortical_2010\')

drive_letter = 'F';

load overall_allcells_coherence_clustered 
load overall_EC_UDS_sep
load corresponding_lfp_state_data_hipp_lfp
load corresponding_lfp_state_data
% good_hipp_lfps(good_hipp_lfps > 109) = [];

%% restrict analysis to good hippocampal LFPs (using best of either LF3 or
%% LF2r)
% bad_lf2r = find(lf2r_kl_lf==0);
% bad_lf3 = find(lf3_kl_lf==0);
% use_lf2r(ismember(use_lf2r,bad_lf2r)) = [];
% use_lf3(ismember(use_lf3,bad_lf3)) = [];
% used_hipp_lfp = union(use_lf2r,use_lf3);
% sess_data = sess_data(used_hipp_lfp);
% 
% use_lf2r = find(ismember(used_hipp_lfp,use_lf2r));
% use_lf3 = find(ismember(used_hipp_lfp,use_lf3));
% 
% lf2r_kl_lf = lf2r_kl_lf(used_hipp_lfp);
% lf3_kl_lf = lf3_kl_lf(used_hipp_lfp);
% lf8_kl_lf = lf8_kl_lf(used_hipp_lfp);
% mp_kl_lf = mp_kl_lf(used_hipp_lfp);
% hipp_kl_lf = zeros(size(used_hipp_lfp));
% hipp_kl_lf(use_lf2r) = lf2r_kl_lf(use_lf2r);
% hipp_kl_lf(use_lf3) = lf3_kl_lf(use_lf3);
% 
% pers_fract_across_cycles_lf2r = pers_fract_across_cycles_lf2r(used_hipp_lfp);
% pers_fract_within_cycle_lf2r = pers_fract_within_cycle_lf2r(used_hipp_lfp);
% pers_fract_across_cycles_lf3 = pers_fract_across_cycles_lf3(used_hipp_lfp);
% pers_fract_within_cycle_lf3 = pers_fract_within_cycle_lf3(used_hipp_lfp);
% pers_fract_across_cycles = pers_fract_across_cycles(used_hipp_lfp(used_hipp_lfp < 116));
% pers_fract_within_cycle = pers_fract_within_cycle(used_hipp_lfp(used_hipp_lfp < 116));
% mean_up_lag_lf2r = mean_up_lag_lf2r(used_hipp_lfp);
% mean_down_lag_lf2r = mean_down_lag_lf2r(used_hipp_lfp);
% mean_up_lag_lf3 = mean_up_lag_lf3(used_hipp_lfp);
% mean_down_lag_lf3 = mean_down_lag_lf3(used_hipp_lfp);
% mean_up_lag = mean_up_lag(used_hipp_lfp(used_hipp_lfp < 116));
% mean_down_lag = mean_down_lag(used_hipp_lfp(used_hipp_lfp < 116));
% median_up_lag_lf2r = median_up_lag_lf2r(used_hipp_lfp);
% median_down_lag_lf2r = median_down_lag_lf2r(used_hipp_lfp);
% median_up_lag_lf3 = median_up_lag_lf3(used_hipp_lfp);
% median_down_lag_lf3 = median_down_lag_lf3(used_hipp_lfp);
% median_up_lag = median_up_lag(used_hipp_lfp(used_hipp_lfp < 116));
% median_down_lag = median_down_lag(used_hipp_lfp(used_hipp_lfp < 116));
% 
% pers_fract_across_cycles_hipp = nan(size(used_hipp_lfp));
% pers_fract_across_cycles_hipp(use_lf2r) = pers_fract_across_cycles_lf2r(use_lf2r);
% pers_fract_across_cycles_hipp(use_lf3) = pers_fract_across_cycles_lf3(use_lf3);
% pers_fract_within_cycle_hipp = nan(size(used_hipp_lfp));
% pers_fract_within_cycle_hipp(use_lf2r) = pers_fract_within_cycle_lf2r(use_lf2r);
% pers_fract_within_cycle_hipp(use_lf3) = pers_fract_within_cycle_lf3(use_lf3);
% mean_up_lag_hipp = nan(size(used_hipp_lfp));
% mean_up_lag_hipp(use_lf2r) = mean_up_lag_lf2r(use_lf2r);
% mean_up_lag_hipp(use_lf3) = mean_up_lag_lf3(use_lf3);
% median_up_lag_hipp = nan(size(used_hipp_lfp));
% median_up_lag_hipp(use_lf2r) = median_up_lag_lf2r(use_lf2r);
% median_up_lag_hipp(use_lf3) = median_up_lag_lf3(use_lf3);
% mean_down_lag_hipp = nan(size(used_hipp_lfp));
% mean_down_lag_hipp(use_lf2r) = mean_down_lag_lf2r(use_lf2r);
% mean_down_lag_hipp(use_lf3) = mean_down_lag_lf3(use_lf3);
% median_down_lag_hipp = nan(size(used_hipp_lfp));
% median_down_lag_hipp(use_lf2r) = median_down_lag_lf2r(use_lf2r);
% median_down_lag_hipp(use_lf3) = median_down_lag_lf3(use_lf3);

%% Use only LF2
bad_lf2r = find(lf2r_kl_lf==0);
good_hipp_lf2r(ismember(good_hipp_lf2r,bad_lf2r)) = [];
sess_data = sess_data(good_hipp_lf2r);

lf2r_kl_lf = lf2r_kl_lf(good_hipp_lf2r);
lf8_kl_lf = lf8_kl_lf(good_hipp_lf2r);
mp_kl_lf = mp_kl_lf(good_hipp_lf2r);
hipp_kl_lf = lf2r_kl_lf;
uds_time_hipp = uds_time_2r(good_hipp_lf2r);
uds_time_lf8 = uds_time_lf8(good_hipp_lf2r);

pers_fract_across_cycles_hipp = pers_fract_across_cycles_lf2r(good_hipp_lf2r);
pers_fract_within_cycle_hipp = pers_fract_within_cycle_lf2r(good_hipp_lf2r);
mean_up_lag_hipp = mean_up_lag_lf2r(good_hipp_lf2r);
mean_down_lag_hipp = mean_down_lag_lf2r(good_hipp_lf2r);
mean_up_lag = mean_up_lag(good_hipp_lf2r(good_hipp_lf2r < 116));
mean_down_lag = mean_down_lag(good_hipp_lf2r(good_hipp_lf2r < 116));
median_up_lag_hipp = median_up_lag_lf2r(good_hipp_lf2r);
median_down_lag_hipp = median_down_lag_lf2r(good_hipp_lf2r);
median_up_lag = median_up_lag(good_hipp_lf2r(good_hipp_lf2r < 116));
median_down_lag = median_down_lag(good_hipp_lf2r(good_hipp_lf2r < 116));

lf2r_updur_lfpc = lf2r_updur_lfpc(good_hipp_lf2r);
mp_updur_lfpc = mp_updur_lfpc(good_hipp_lf2r < 116);

%%
mec = find_struct_field_vals(sess_data,'region','MEC');
layer3 = find_struct_field_vals(sess_data,'layer','3');
layer2 = find_struct_field_vals(sess_data,'layer','2');
layer23 = find_struct_field_vals(sess_data,'layer','23');
lec = find_struct_field_vals(sess_data,'region','LEC');
l3mec = intersect(mec,layer3);
l2mec = intersect(mec,layer2);
l3lec = intersect(lec,layer3);
l3mec(24:end) = [];
l23mec = intersect(mec,layer23);
l23mec = unique([l2mec l3mec l23mec]);

hipp_sep_threshold = prctile(lf8_kl_lf,10);
%% plot LF2 separability vs LF8 separability
%(ignoring 0 points which arise when there aren't sufficient length of data
%with UDS
figure
plot(lf8_kl_lf,hipp_kl_lf,'.')
line([0 35],[0 35],'color','k')
xl = xlim();
line(xl,[hipp_sep_threshold hipp_sep_threshold],'color','k','linestyle','--')
xlabel('LF8 UDS separability','fontsize',14)
ylabel('LF2r UDS separability','fontsize',14)

figure
plot(uds_time_lf8*5,uds_time_hipp,'.')
line([0 1800],[0 1800],'color','k')
xl = xlim();
xlabel('LF8 UDS duration (s)','fontsize',14)
ylabel('LF2r UDS duration (s)','fontsize',14)

%% plot level of LF2 persistent activity vs level of MP persistent activity
figure
plot(pers_fract_across_cycles_hipp(l3mec),pers_fract_across_cycles(l3mec),'.')
hold on
usedmec = l3mec(hipp_kl_lf(l3mec) > hipp_sep_threshold);
plot(pers_fract_across_cycles_hipp(usedmec),pers_fract_across_cycles(usedmec),'o')
plot(pers_fract_across_cycles_hipp(l3lec),pers_fract_across_cycles(l3lec),'r.')
usedlec = l3lec(hipp_kl_lf(l3lec) > hipp_sep_threshold);
plot(pers_fract_across_cycles_hipp(usedlec),pers_fract_across_cycles(usedlec),'ro')
line([0 0.6],[0 0.6],'color','k')
xlim([0 0.4])
xlabel('Percent type 2 persistence (hippocampal LFP)','fontsize',14)
ylabel('Percent type 2 persistence (MP)','fontsize',14)
legend('L3MEC','L3MEC clear Hipp UDS','L3LEC','L3LEC clear Hipp UDS')

%% plot avg UP Lag
figure
plot(median_up_lag_hipp(l3mec),median_up_lag(l3mec),'.')
hold on
plot(median_up_lag_hipp(usedmec),median_up_lag(usedmec),'o')
plot(median_up_lag_hipp(l3lec),median_up_lag(l3lec),'r.')
plot(median_up_lag_hipp(usedlec),median_up_lag(usedlec),'ro')
line([-0.1 0.6],[-0.1 0.6],'color','k')
xlabel('Median up lag (s) (hippocampal LFP)','fontsize',14)
ylabel('Median up lag (s) (MP)','fontsize',14)
legend('L3MEC','L3MEC clear Hipp UDS','L3LEC','L3LEC clear Hipp UDS')

%% plot avg DOWN Lag
figure
plot(median_down_lag_hipp(l3mec),median_down_lag(l3mec),'.')
hold on
plot(median_down_lag_hipp(usedmec),median_down_lag(usedmec),'o')
plot(median_down_lag_hipp(l3lec),median_down_lag(l3lec),'r.')
plot(median_down_lag_hipp(usedlec),median_down_lag(usedlec),'ro')
line([0 1.2],[0 1.2],'color','k')
xlabel('Median down lag (s) (hippocampal LFP)','fontsize',14)
ylabel('Median down lag (s) (MP)','fontsize',14)
legend('L3MEC','L3MEC clear Hipp UDS','L3LEC','L3LEC clear Hipp UDS')


%%
for i = 1:length(mp_updur_lfpc)
   mean_updur_lpc(i) = nanmean(mp_updur_lfpc{i});
   mean_updur_lpc_lf2r(i) = nanmean(lf2r_updur_lfpc{i});
end

figure
plot(mean_updur_lpc_lf2r(l3mec),mean_updur_lpc(l3mec),'.')
hold on
plot(mean_updur_lpc_lf2r(usedmec),mean_updur_lpc(usedmec),'o')
plot(mean_updur_lpc_lf2r(l3lec),mean_updur_lpc(l3lec),'r.')
plot(mean_updur_lpc_lf2r(usedlec),mean_updur_lpc(usedlec),'ro')



