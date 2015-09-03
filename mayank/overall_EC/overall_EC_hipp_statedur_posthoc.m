clear all
% close all
%%
load F:\WC_Germany\overall_EC\overall_allcells_dir
addpath('F:\Code\WC_anal\general\')
addpath('F:\WC_Germany\Overall_EC\')
addpath('F:\Code\Chronux\spectral_analysis\continuous\')
addpath('F:\WC_Germany\hsmm_state_detection\\')
addpath('F:\WC_Germany\parietal_cortical_2010\')

drive_letter = 'F';

cd F:\WC_Germany\overall_EC
load ./overall_allcells_coherence_clustered 
load ./overall_EC_UDS_sep
load ./corresponding_lfp_state_data_hipp_lfp
load ./corresponding_lfp_state_data
load ./state_dur_stats
% good_hipp_lfps(good_hipp_lfps > 109) = [];

%% Use only LF2
bad_lf2r = find(lf2r_kl_lf==0);
good_hipp_lf2r(ismember(good_hipp_lf2r,bad_lf2r)) = [];

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
clear_hipp_lf2r = find(lf2r_kl_lf(good_hipp_lf2r) > hipp_sep_threshold);

%% For up duration distributions
maxdur = 10;
up_range = [0.1 maxdur];
down_range = [0.1 maxdur];
numBins = 60;

for d = 1:length(l3mec)
    [l3mec_up_hist(d,:),up_grid] = log_hist(mp_state_durations{l3mec(d)}{2}(mp_state_durations{l3mec(d)}{2} <= maxdur),up_range,numBins);
end
for d = 1:length(l3lec)
    [l3lec_up_hist(d,:),up_grid] = log_hist(mp_state_durations{l3lec(d)}{2}(mp_state_durations{l3lec(d)}{2} <= maxdur),up_range,numBins);
end
for d = 1:length(good_hipp_lf2r)
    [hipp_up_hist(d,:),up_grid] = log_hist(lf2r_state_durations{good_hipp_lf2r(d)}{2}(lf2r_state_durations{good_hipp_lf2r(d)}{2} <= maxdur),up_range,numBins);
end

figure
h=errorbar(up_grid,nanmean(l3mec_up_hist),nanstd(l3mec_up_hist)/sqrt(length(l3mec)));
errorbar_tick(h,.01,'units');
hold on
h=errorbar(up_grid,nanmean(l3lec_up_hist),nanstd(l3lec_up_hist)/sqrt(length(l3lec)),'g');
errorbar_tick(h,.01,'units');
% h=errorbar(up_grid,nanmean(hipp_up_hist),nanstd(hipp_up_hist)/sqrt(length(good_hipp_lf2r)),'k');
% errorbar_tick(h,.01,'units');
h=errorbar(up_grid,nanmean(hipp_up_hist(clear_hipp_lf2r,:)),nanstd(hipp_up_hist(clear_hipp_lf2r,:))/sqrt(length(clear_hipp_lf2r)),'k');
errorbar_tick(h,.01,'units');
set(gca,'yscale','log')
xlim([0 maxdur])
ylim([1e-4 0.15])
xlabel('Up state duration (s)','fontsize',14)
ylabel('Relative Frequency','fontsize',14)

%% For down duration distributions
maxdur = 10;
up_range = [0.1 maxdur];
down_range = [0.1 maxdur];
numBins = 60;

for d = 1:length(l3mec)
    [l3mec_up_hist(d,:),up_grid] = log_hist(mp_state_durations{l3mec(d)}{1}(mp_state_durations{l3mec(d)}{1} <= maxdur),up_range,numBins);
end
for d = 1:length(l3lec)
    [l3lec_up_hist(d,:),up_grid] = log_hist(mp_state_durations{l3lec(d)}{1}(mp_state_durations{l3lec(d)}{1} <= maxdur),up_range,numBins);
end
for d = 1:length(good_hipp_lf2r)
    [hipp_up_hist(d,:),up_grid] = log_hist(lf2r_state_durations{good_hipp_lf2r(d)}{1}(lf2r_state_durations{good_hipp_lf2r(d)}{1} <= maxdur),up_range,numBins);
end

figure
h=errorbar(up_grid,nanmean(l3mec_up_hist),nanstd(l3mec_up_hist)/sqrt(length(l3mec)));
errorbar_tick(h,.01,'units');
hold on
h=errorbar(up_grid,nanmean(l3lec_up_hist),nanstd(l3lec_up_hist)/sqrt(length(l3lec)),'g');
errorbar_tick(h,.01,'units');
% h=errorbar(up_grid,nanmean(hipp_up_hist),nanstd(hipp_up_hist)/sqrt(length(good_hipp_lf2r)),'k');
% errorbar_tick(h,.01,'units');
h=errorbar(up_grid,nanmean(hipp_up_hist(clear_hipp_lf2r,:)),nanstd(hipp_up_hist(clear_hipp_lf2r,:))/sqrt(length(clear_hipp_lf2r)),'k');
errorbar_tick(h,.01,'units');
set(gca,'yscale','log')
xlim([0 maxdur])
ylim([1e-4 0.15])
xlabel('Down state duration (s)','fontsize',14)
ylabel('Relative Frequency','fontsize',14)


%%
hist_range = linspace(0,9,1000);
binsize = hist_range(2)-hist_range(1);
for i = 1:length(l3mec)
    l3mec_hist(i,:) = hist(mp_updur_lfpc{l3mec(i)},hist_range);
    l3mec_hist(i,:) = l3mec_hist(i,:)/sum(l3mec_hist(i,:))/binsize;
    l3mec_hist(i,:) = jmm_smooth_1d_cor(l3mec_hist(i,:),20);
end
for i = 1:length(l3lec)
    l3lec_hist(i,:) = hist(mp_updur_lfpc{l3lec(i)},hist_range);
    l3lec_hist(i,:) = l3lec_hist(i,:)/sum(l3lec_hist(i,:))/binsize;
    l3lec_hist(i,:) = jmm_smooth_1d_cor(l3lec_hist(i,:),20);
end
for i = 1:length(good_hipp_lf2r)
    hipp_hist(i,:) = hist(lf2r_updur_lfpc{good_hipp_lf2r(i)},hist_range);
    hipp_hist(i,:) = hipp_hist(i,:)/sum(hipp_hist(i,:))/binsize;
    hipp_hist(i,:) = jmm_smooth_1d_cor(hipp_hist(i,:),20);
end


figure
h=errorbar(hist_range,nanmean(l3mec_hist),nanstd(l3mec_hist)/sqrt(length(l3mec)));
errorbar_tick(h,.01,'units');
hold on
h=errorbar(hist_range,nanmean(l3lec_hist),nanstd(l3lec_hist)/sqrt(length(l3lec)),'g');
errorbar_tick(h,.01,'units');
h=errorbar(hist_range,nanmean(hipp_hist(clear_hipp_lf2r,:)),nanstd(hipp_hist(clear_hipp_lf2r,:))/sqrt(length(clear_hipp_lf2r)),'k');
errorbar_tick(h,.01,'units');
set(gca,'yscale','log')
xlim([0 8])
ylim([1e-4 3])
xlabel('Up state duration (LFP Cycles)','fontsize',14)
ylabel('Probability density','fontsize',14)

%% UP transition lag distributions
lag_range = linspace(-3,3,1000);
dlag = lag_range(2)-lag_range(1);
for d = 1:length(l3mec)
    l3mec_hist(d,:) = histc(lfp_up_lag{l3mec(d)},lag_range);
    l3mec_hist(d,:) = l3mec_hist(d,:)/sum(l3mec_hist(d,:))/dlag;
    sm_l3mec_hist(d,:) = jmm_smooth_1d_cor(l3mec_hist(d,:),10);
end
for d = 1:length(l3lec)
    l3lec_hist(d,:) = histc(lfp_up_lag{l3lec(d)},lag_range);
    l3lec_hist(d,:) = l3lec_hist(d,:)/sum(l3lec_hist(d,:))/dlag;
    sm_l3lec_hist(d,:) = jmm_smooth_1d_cor(l3lec_hist(d,:),10);
end
for d = 1:length(good_hipp_lf2r)
    hipp_hist(d,:) = histc(lf2r_up_lag{good_hipp_lf2r(d)},lag_range);
    hipp_hist(d,:) = hipp_hist(d,:)/sum(hipp_hist(d,:))/dlag;
    sm_hipp_hist(d,:) = jmm_smooth_1d_cor(hipp_hist(d,:),10);
end

figure
h=errorbar(lag_range,nanmean(sm_l3mec_hist),nanstd(sm_l3mec_hist)/sqrt(length(l3mec)));
errorbar_tick(h,.01,'units');
hold on
h=errorbar(lag_range,nanmean(sm_l3lec_hist),nanstd(sm_l3lec_hist)/sqrt(length(l3lec)),'g');
errorbar_tick(h,.01,'units');
h=errorbar(lag_range,nanmean(sm_hipp_hist(clear_hipp_lf2r,:)),nanstd(sm_hipp_hist(clear_hipp_lf2r,:))/sqrt(length(clear_hipp_lf2r)),'k');
errorbar_tick(h,.01,'units');
xlim([-0.5 1.5])
% ylim([0 3])
line([0 0],[0 3],'Color','k')
xlabel('Up-Transition Lag (s)','fontsize',14)
ylabel('Probability','fontsize',14)
legend('L3MEC','L3LEC','Hipp LFP')

%% DOWN transition lag distributions
lag_range = linspace(-1,5,1000);
dlag = lag_range(2)-lag_range(1);
for d = 1:length(l3mec)
    l3mec_hist(d,:) = histc(lfp_down_lag{l3mec(d)},lag_range);
    l3mec_hist(d,:) = l3mec_hist(d,:)/sum(l3mec_hist(d,:))/dlag;
    sm_l3mec_hist(d,:) = jmm_smooth_1d_cor(l3mec_hist(d,:),10);
end
for d = 1:length(l3lec)
    l3lec_hist(d,:) = histc(lfp_down_lag{l3lec(d)},lag_range);
    l3lec_hist(d,:) = l3lec_hist(d,:)/sum(l3lec_hist(d,:))/dlag;
    sm_l3lec_hist(d,:) = jmm_smooth_1d_cor(l3lec_hist(d,:),10);
end
for d = 1:length(good_hipp_lf2r)
    hipp_hist(d,:) = histc(lf2r_down_lag{good_hipp_lf2r(d)},lag_range);
    hipp_hist(d,:) = hipp_hist(d,:)/sum(hipp_hist(d,:))/dlag;
    sm_hipp_hist(d,:) = jmm_smooth_1d_cor(hipp_hist(d,:),10);
end

figure
h=errorbar(lag_range,nanmean(sm_l3mec_hist),nanstd(sm_l3mec_hist)/sqrt(length(l3mec)));
errorbar_tick(h,.01,'units');
hold on
h=errorbar(lag_range,nanmean(sm_l3lec_hist),nanstd(sm_l3lec_hist)/sqrt(length(l3lec)),'g');
errorbar_tick(h,.01,'units');
h=errorbar(lag_range,nanmean(sm_hipp_hist(clear_hipp_lf2r,:)),nanstd(sm_hipp_hist(clear_hipp_lf2r,:))/sqrt(length(clear_hipp_lf2r)),'k');
errorbar_tick(h,.01,'units');
xlim([-0.5 2.5])
% ylim([0 3])
line([0 0],[0 2],'Color','k')
xlabel('Down-Transition Lag (s)','fontsize',14)
ylabel('Probability','fontsize',14)
legend('L3MEC','L3LEC','Hipp LFP')
