clear all
close all
%%
load E:\WC_Germany\overall_EC\overall_allcells_dir
addpath('E:\Code\WC_anal\general\')
addpath('E:\WC_Germany\Overall_EC\')
addpath('E:\Code\Chronux\spectral_analysis\continuous\')
addpath('E:\WC_Germany\hsmm_state_detection\\')
addpath('E:\WC_Germany\parietal_cortical_2010\')

drive_letter = 'E';
dsf = 16;
Fsd = 2016/dsf;

load overall_EC_trig_lf2lfhf_layer3_data
load overall_allcells_coherence_clustered 
good_hipp_lf2r(good_hipp_lf2r > 109) = [];

%% restrict analysis to good hippocampal LFPs
sess_data = sess_data(good_hipp_lf2r); 
mp_utrig_mp = mp_utrig_mp(good_hipp_lf2r,:);
mp_utrig_lf8 = mp_utrig_lf8(good_hipp_lf2r,:);
mp_utrig_lf2rlf = mp_utrig_lf2rlf(good_hipp_lf2r,:);
mp_utrig_lf2rhf = mp_utrig_lf2rhf(good_hipp_lf2r,:);
mp_utrig_lf2hf = mp_utrig_lf2hf(good_hipp_lf2r,:);
mp_utrig_spk = mp_utrig_spk(good_hipp_lf2r,:)*Fsd;
mp_dtrig_mp = mp_dtrig_mp(good_hipp_lf2r,:);
mp_dtrig_lf8 = mp_dtrig_lf8(good_hipp_lf2r,:);
mp_dtrig_lf2rlf = mp_dtrig_lf2rlf(good_hipp_lf2r,:);
mp_dtrig_lf2rhf = mp_dtrig_lf2rhf(good_hipp_lf2r,:);
mp_dtrig_lf2hf = mp_dtrig_lf2hf(good_hipp_lf2r,:);
mp_dtrig_spk = mp_dtrig_spk(good_hipp_lf2r,:)*Fsd;
lf8_utrig_lf8 = lf8_utrig_lf8(good_hipp_lf2r,:);
lf8_utrig_lf2rlf = lf8_utrig_lf2rlf(good_hipp_lf2r,:);
lf8_utrig_lf2rhf = lf8_utrig_lf2rhf(good_hipp_lf2r,:);
lf8_utrig_lf2hf = lf8_utrig_lf2hf(good_hipp_lf2r,:);
lf8_utrig_mp = lf8_utrig_mp(good_hipp_lf2r,:);
lf8_utrig_spk = lf8_utrig_spk(good_hipp_lf2r,:)*Fsd;
lf8_dtrig_lf8 = lf8_dtrig_lf8(good_hipp_lf2r,:);
lf8_dtrig_lf2rlf = lf8_dtrig_lf2rlf(good_hipp_lf2r,:);
lf8_dtrig_lf2rhf = lf8_dtrig_lf2rhf(good_hipp_lf2r,:);
lf8_dtrig_lf2hf = lf8_dtrig_lf2hf(good_hipp_lf2r,:);
lf8_dtrig_spk = lf8_dtrig_spk(good_hipp_lf2r,:)*Fsd;
lf8_dtrig_mp = lf8_dtrig_mp(good_hipp_lf2r,:);
% lf8_utrig_lf2r_cond_mdown = lf8_utrig_lf2r_cond_mdown(good_hipp_lf2r,:);
% lf8_utrig_spk_cond_mdown = lf8_utrig_spk_cond_mdown(good_hipp_lf2r,:)*Fsd;
% lf8_utrig_mp_cond_mdown = lf8_utrig_mp_cond_mdown(good_hipp_lf2r,:);
% lf8_dtrig_lf2r_cond_mup = lf8_dtrig_lf2r_cond_mup(good_hipp_lf2r,:);
% lf8_dtrig_spk_cond_mup = lf8_dtrig_spk_cond_mup(good_hipp_lf2r,:)*Fsd;
% lf8_dtrig_mp_cond_mup = lf8_dtrig_mp_cond_mup(good_hipp_lf2r,:);
% lf8_utrig_lf2r_cond_mup = lf8_utrig_lf2r_cond_mup(good_hipp_lf2r,:);
% lf8_utrig_spk_cond_mup = lf8_utrig_spk_cond_mup(good_hipp_lf2r,:)*Fsd;
% lf8_utrig_mp_cond_mup = lf8_utrig_mp_cond_mup(good_hipp_lf2r,:);
% 
%%
mec = find_struct_field_vals(sess_data,'region','MEC');
layer3 = find_struct_field_vals(sess_data,'layer','3');
layer2 = find_struct_field_vals(sess_data,'layer','2');
layer23 = find_struct_field_vals(sess_data,'layer','23');
layer5 = find_struct_field_vals(sess_data,'layer','56');
lec = find_struct_field_vals(sess_data,'region','LEC');
l3mec = intersect(mec,layer3);
l2mec = intersect(mec,layer2);
l3lec = intersect(lec,layer3);
l3mec(24:end) = [];
l5mec = intersect(mec,layer5);
l23mec = intersect(mec,layer23);
l23mec = unique([l2mec l3mec l23mec]);
ca1 = find_struct_field_vals(sess_data,'region','ca1');
ca3 = find_struct_field_vals(sess_data,'region','ca3');
dg = find_struct_field_vals(sess_data,'region','dg');
pyr = find_struct_field_vals(sess_data,'cell_type','pyr');
ca1p = intersect(ca1,pyr);
ca1i = setdiff(ca1,pyr);

%
brown = [0.3 0.2 0.3];

%% CREATE MP UP-TRIG AVG COMPARISON
figure
h = errorbar(lags/Fsd,nanmean(mp_utrig_mp(l3mec,:)),nanstd(mp_utrig_mp(l3mec,:))/sqrt(length(l3mec)));
errorbar_tick(h,.01,'units');
hold on
h = errorbar(lags/Fsd,nanmean(mp_utrig_lf8(l3mec,:)),nanstd(mp_utrig_lf8(l3mec,:))/sqrt(length(l3mec)),'r');
errorbar_tick(h,.01,'units');
h = errorbar(lags/Fsd,nanmean(mp_utrig_lf2rlf(l3mec,:)),nanstd(mp_utrig_lf2rlf(l3mec,:))/sqrt(length(l3mec)),'g');
errorbar_tick(h,.01,'units');
h = errorbar(lags/Fsd,nanmean(mp_utrig_lf2rhf(l3mec,:)),nanstd(mp_utrig_lf2rhf(l3mec,:))/sqrt(length(l3mec)),'y');
errorbar_tick(h,.01,'units');
h = errorbar(lags/Fsd,nanmean(mp_utrig_lf2hf(l3mec,:)),nanstd(mp_utrig_lf2hf(l3mec,:))/sqrt(length(l3mec)),'k');
errorbar_tick(h,.01,'units');
h = errorbar(lags/Fsd,nanmean(mp_utrig_spk(l3mec,:))/8,nanstd(mp_utrig_spk(l3mec,:))/8/sqrt(length(l3mec)),'c');
errorbar_tick(h,.01,'units');
xlim([-0.5 1]), grid
title('L3MEC Up Transition')
legend('MP','LF8','LF2r-LF','LF2r-HF','LF2-HF','Firing rate')
yl = ylim();
line([0 0],yl,'color','k')
xlabel('Time lag (s)','fontsize',14)
ylabel('Amplitude (z)','fontsize',14)


figure
h = errorbar(lags/Fsd,nanmean(mp_utrig_lf2r(l3mec,:)),nanstd(mp_utrig_lf2r(l3mec,:))/sqrt(length(l3mec)));
errorbar_tick(h,.01,'units');
hold on
h = errorbar(lags/Fsd,nanmean(mp_utrig_lf8(l3mec,:)),nanstd(mp_utrig_lf8(l3mec,:))/sqrt(length(l3mec)),'r');
errorbar_tick(h,.01,'units');
h = errorbar(lags/Fsd,nanmean(mp_utrig_lf2r(l3lec,:)),nanstd(mp_utrig_lf2r(l3lec,:))/sqrt(length(l3lec)),'k');
errorbar_tick(h,.01,'units');
h = errorbar(lags/Fsd,nanmean(mp_utrig_lf8(l3lec,:)),nanstd(mp_utrig_lf8(l3lec,:))/sqrt(length(l3lec)),'c');
errorbar_tick(h,.01,'units');
xlim([-0.5 1]), grid
legend('MEC U-trig LF2r','MEC U-trig LF8','LEC U-trig LF2r','LEC U-trig LF8')
yl = ylim();
line([0 0],yl,'color','k')
xlabel('Time lag (s)','fontsize',14)
ylabel('Amplitude (z)','fontsize',14)

figure
h = errorbar(lags/Fsd,nanmean(mp_utrig_mp(l3lec,:)),nanstd(mp_utrig_mp(l3lec,:))/sqrt(length(l3lec)));
errorbar_tick(h,.01,'units');
hold on
h = errorbar(lags/Fsd,nanmean(mp_utrig_lf8(l3lec,:)),nanstd(mp_utrig_lf8(l3lec,:))/sqrt(length(l3lec)),'r');
errorbar_tick(h,.01,'units');
h = errorbar(lags/Fsd,nanmean(mp_utrig_lf2rlf(l3lec,:)),nanstd(mp_utrig_lf2rlf(l3lec,:))/sqrt(length(l3lec)),'g');
errorbar_tick(h,.01,'units');
h = errorbar(lags/Fsd,nanmean(mp_utrig_lf2rhf(l3lec,:)),nanstd(mp_utrig_lf2rhf(l3lec,:))/sqrt(length(l3lec)),'y');
errorbar_tick(h,.01,'units');
h = errorbar(lags/Fsd,nanmean(mp_utrig_lf2hf(l3lec,:)),nanstd(mp_utrig_lf2hf(l3lec,:))/sqrt(length(l3lec)),'k');
errorbar_tick(h,.01,'units');
h = errorbar(lags/Fsd,nanmean(mp_utrig_spk(l3lec,:))/8,nanstd(mp_utrig_spk(l3lec,:))/8/sqrt(length(l3lec)),'c');
errorbar_tick(h,.01,'units');
xlim([-0.5 1]), grid
title('L3LEC Up Transition')
legend('MP','LF8','LF2r-LF','LF2r-HF','LF2-HF','Firing rate')
yl = ylim();
line([0 0],yl,'color','k')
xlabel('Time lag (s)','fontsize',14)
ylabel('Amplitude (z)','fontsize',14)

%% CREATE LFP UP-TRIG AVG COMPARISON
figure
h = errorbar(lags/Fsd,nanmean(lf8_utrig_mp(l3mec,:)),nanstd(lf8_utrig_mp(l3mec,:))/sqrt(length(l3mec)));
errorbar_tick(h,.01,'units');
hold on
h = errorbar(lags/Fsd,nanmean(lf8_utrig_mp(l3lec,:)),nanstd(lf8_utrig_mp(l3lec,:))/sqrt(length(l3lec)),'g');
errorbar_tick(h,.01,'units');
h = errorbar(lags/Fsd,nanmean(lf8_utrig_lf8),nanstd(lf8_utrig_lf8)/sqrt(length(sess_data)),'r');
errorbar_tick(h,.01,'units');
h = errorbar(lags/Fsd,nanmean(lf8_utrig_lf2rlf),nanstd(lf8_utrig_lf2rlf)/sqrt(length(sess_data)),'k');
errorbar_tick(h,.01,'units');
xlim([-0.5 1]), grid
title('LF8 Up Transition')
legend('L3MEC','L3LEC','LF8','LF2r')
yl = ylim();
line([0 0],yl,'color','k')
xlabel('Time lag (s)','fontsize',14)
ylabel('Amplitude (z)','fontsize',14)


figure
h = errorbar(lags/Fsd,nanmean(lf8_utrig_mp(l3mec,:)),nanstd(lf8_utrig_mp(l3mec,:))/sqrt(length(l3mec)));
errorbar_tick(h,.01,'units');
hold on
h=errorbar(lags/Fsd,nanmean(lf8_utrig_mp(l2mec,:)),nanstd(lf8_utrig_mp(l2mec,:))/sqrt(length(l2mec)),'color',[0.7 0.1 0.2]);
errorbar_tick(h,.01,'units')
h=errorbar(lags/Fsd,nanmean(lf8_utrig_mp(l5mec,:)),nanstd(lf8_utrig_mp(l5mec,:))/sqrt(length(l5mec)),'color',[0.2 0.5 0.4]);
errorbar_tick(h,.01,'units')
h=errorbar(lags/Fsd,nanmean(lf8_utrig_mp(l3lec,:)),nanstd(lf8_utrig_mp(l3lec,:))/sqrt(length(l3lec)),'color',[0.8 0.1 0.1]);
errorbar_tick(h,.01,'units')
h=errorbar(lags/Fsd,nanmean(lf8_utrig_mp(ca1p,:)),nanstd(lf8_utrig_mp(ca1p,:))/sqrt(length(ca1p)),'k');
errorbar_tick(h,.01,'units')
h=errorbar(lags/Fsd,nanmean(lf8_utrig_mp(ca1i,:)),nanstd(lf8_utrig_mp(ca1i,:))/sqrt(length(ca1i)),'g');
errorbar_tick(h,.01,'units')
h=errorbar(lags/Fsd,nanmean(lf8_utrig_mp(dg,:)),nanstd(lf8_utrig_mp(dg,:))/sqrt(length(dg)),'c');
errorbar_tick(h,.01,'units')
h=errorbar(lags/Fsd,nanmean(lf8_utrig_mp(ca3,:)),nanstd(lf8_utrig_mp(ca3,:))/sqrt(length(ca3)),'y');
errorbar_tick(h,.01,'units')
h=errorbar(lags/Fsd,nanmean(lf8_utrig_lf8),nanstd(lf8_utrig_lf8)/sqrt(length(sess_data)),'r');
errorbar_tick(h,.01,'units')
xlim([-0.5 1]), grid
title('LF8 Up Transition')
legend('L3MEC','L2MEC','L5MEC','L3LEC','CA1p','CA1i','DG','CA3','LF8')
yl = ylim();
line([0 0],yl,'color','k')
xlabel('Time lag (s)','fontsize',14)
ylabel('Amplitude (z)','fontsize',14)

figure
h = errorbar(lags/Fsd,nanmean(lf8_dtrig_mp(l3mec,:)),nanstd(lf8_dtrig_mp(l3mec,:))/sqrt(length(l3mec)));
errorbar_tick(h,.01,'units');
hold on
h=errorbar(lags/Fsd,nanmean(lf8_dtrig_mp(l2mec,:)),nanstd(lf8_dtrig_mp(l2mec,:))/sqrt(length(l2mec)),'color',[0.7 0.1 0.2]);
errorbar_tick(h,.01,'units')
h=errorbar(lags/Fsd,nanmean(lf8_dtrig_mp(l5mec,:)),nanstd(lf8_dtrig_mp(l5mec,:))/sqrt(length(l5mec)),'color',[0.2 0.5 0.4]);
errorbar_tick(h,.01,'units')
h=errorbar(lags/Fsd,nanmean(lf8_dtrig_mp(l3lec,:)),nanstd(lf8_dtrig_mp(l3lec,:))/sqrt(length(l3lec)),'color',[0.8 0.1 0.1]);
errorbar_tick(h,.01,'units')
h=errorbar(lags/Fsd,nanmean(lf8_dtrig_mp(ca1p,:)),nanstd(lf8_dtrig_mp(ca1p,:))/sqrt(length(ca1p)),'k');
errorbar_tick(h,.01,'units')
h=errorbar(lags/Fsd,nanmean(lf8_dtrig_mp(ca1i,:)),nanstd(lf8_dtrig_mp(ca1i,:))/sqrt(length(ca1i)),'g');
errorbar_tick(h,.01,'units')
h=errorbar(lags/Fsd,nanmean(lf8_dtrig_mp(dg,:)),nanstd(lf8_dtrig_mp(dg,:))/sqrt(length(dg)),'c');
errorbar_tick(h,.01,'units')
h=errorbar(lags/Fsd,nanmean(lf8_dtrig_mp(ca3,:)),nanstd(lf8_dtrig_mp(ca3,:))/sqrt(length(ca3)),'y');
errorbar_tick(h,.01,'units')
h=errorbar(lags/Fsd,nanmean(lf8_dtrig_lf8),nanstd(lf8_dtrig_lf8)/sqrt(length(sess_data)),'r');
errorbar_tick(h,.01,'units')
xlim([-0.5 1]), grid
title('LF8 Down Transition')
legend('L3MEC','L2MEC','L5MEC','L3LEC','CA1p','CA1i','DG','CA3','LF8')
yl = ylim();
line([0 0],yl,'color','k')
xlabel('Time lag (s)','fontsize',14)
ylabel('Amplitude (z)','fontsize',14)

%% CREATE MP DOWN-TRIG AVG COMPARISON
figure
h = errorbar(lags/Fsd,nanmean(mp_dtrig_mp(l3mec,:)),nanstd(mp_dtrig_mp(l3mec,:))/sqrt(length(l3mec)));
errorbar_tick(h,.01,'units');
hold on
h = errorbar(lags/Fsd,nanmean(mp_dtrig_lf8(l3mec,:)),nanstd(mp_dtrig_lf8(l3mec,:))/sqrt(length(l3mec)),'r');
errorbar_tick(h,.01,'units');
h = errorbar(lags/Fsd,nanmean(mp_dtrig_lf2rlf(l3mec,:)),nanstd(mp_dtrig_lf2rlf(l3mec,:))/sqrt(length(l3mec)),'g');
errorbar_tick(h,.01,'units');
h = errorbar(lags/Fsd,nanmean(mp_dtrig_lf2rhf(l3mec,:)),nanstd(mp_dtrig_lf2rhf(l3mec,:))/sqrt(length(l3mec)),'y');
errorbar_tick(h,.01,'units');
h = errorbar(lags/Fsd,nanmean(mp_dtrig_lf2hf(l3mec,:)),nanstd(mp_dtrig_lf2hf(l3mec,:))/sqrt(length(l3mec)),'k');
errorbar_tick(h,.01,'units');
h = errorbar(lags/Fsd,nanmean(mp_dtrig_spk(l3mec,:))/8,nanstd(mp_dtrig_spk(l3mec,:))/8/sqrt(length(l3mec)),'c');
errorbar_tick(h,.01,'units');
xlim([-1 1]), grid
title('L3MEC Up Transition')
legend('MP','LF8','LF2r-LF','LF2r-HF','LF2-HF','Firing rate')
yl = ylim();
line([0 0],yl,'color','k')
xlabel('Time lag (s)','fontsize',14)
ylabel('Amplitude (z)','fontsize',14)

figure
h = errorbar(lags/Fsd,nanmean(mp_dtrig_lf2rlf(l3mec,:)),nanstd(mp_dtrig_lf2rlf(l3mec,:))/sqrt(length(l3mec)));
errorbar_tick(h,.01,'units');
hold on
h = errorbar(lags/Fsd,nanmean(mp_dtrig_lf8(l3mec,:)),nanstd(mp_dtrig_lf8(l3mec,:))/sqrt(length(l3mec)),'r');
errorbar_tick(h,.01,'units');
h = errorbar(lags/Fsd,nanmean(mp_dtrig_lf2rlf(l3lec,:)),nanstd(mp_dtrig_lf2rlf(l3lec,:))/sqrt(length(l3lec)),'k');
errorbar_tick(h,.01,'units');
h = errorbar(lags/Fsd,nanmean(mp_dtrig_lf8(l3lec,:)),nanstd(mp_dtrig_lf8(l3lec,:))/sqrt(length(l3lec)),'c');
errorbar_tick(h,.01,'units');
xlim([-1 1]), grid
legend('MEC D-trig LF2r','MEC D-trig LF8','LEC D-trig LF2r','LEC D-trig LF8')
yl = ylim();
line([0 0],yl,'color','k')
xlabel('Time lag (s)','fontsize',14)
ylabel('Amplitude (z)','fontsize',14)

figure
h = errorbar(lags/Fsd,nanmean(mp_dtrig_mp(l3lec,:)),nanstd(mp_dtrig_mp(l3lec,:))/sqrt(length(l3lec)));
errorbar_tick(h,.01,'units');
hold on
h = errorbar(lags/Fsd,nanmean(mp_dtrig_lf8(l3lec,:)),nanstd(mp_dtrig_lf8(l3lec,:))/sqrt(length(l3lec)),'r');
errorbar_tick(h,.01,'units');
h = errorbar(lags/Fsd,nanmean(mp_dtrig_lf2rlf(l3lec,:)),nanstd(mp_dtrig_lf2rlf(l3lec,:))/sqrt(length(l3lec)),'g');
errorbar_tick(h,.01,'units');
h = errorbar(lags/Fsd,nanmean(mp_dtrig_lf2rhf(l3lec,:)),nanstd(mp_dtrig_lf2rhf(l3lec,:))/sqrt(length(l3lec)),'y');
errorbar_tick(h,.01,'units');
h = errorbar(lags/Fsd,nanmean(mp_dtrig_lf2hf(l3lec,:)),nanstd(mp_dtrig_lf2hf(l3lec,:))/sqrt(length(l3lec)),'k');
errorbar_tick(h,.01,'units');
h = errorbar(lags/Fsd,nanmean(mp_dtrig_spk(l3lec,:))/2,nanstd(mp_dtrig_spk(l3lec,:))/2/sqrt(length(l3lec)),'c');
errorbar_tick(h,.01,'units');
xlim([-1 1]), grid
title('L3LEC Up Transition')
legend('MP','LF8','LF2r-LF','LF2r-HF','LF2-HF','Firing rate')
yl = ylim();
line([0 0],yl,'color','k')
xlabel('Time lag (s)','fontsize',14)
ylabel('Amplitude (z)','fontsize',14)

%% CREATE LFP DOWN-TRIG AVG COMPARISON
figure
h = errorbar(lags/Fsd,nanmean(lf8_dtrig_mp(l3mec,:)),nanstd(lf8_dtrig_mp(l3mec,:))/sqrt(length(l3mec)));
errorbar_tick(h,.01,'units');
hold on
h = errorbar(lags/Fsd,nanmean(lf8_dtrig_mp(l3lec,:)),nanstd(lf8_dtrig_mp(l3lec,:))/sqrt(length(l3lec)),'g');
errorbar_tick(h,.01,'units');
h = errorbar(lags/Fsd,nanmean(lf8_dtrig_lf8),nanstd(lf8_dtrig_lf8)/sqrt(length(sess_data)),'r');
errorbar_tick(h,.01,'units');
h = errorbar(lags/Fsd,nanmean(lf8_dtrig_lf2r),nanstd(lf8_dtrig_lf2r)/sqrt(length(sess_data)),'k');
errorbar_tick(h,.01,'units');
xlim([-0.5 1]), grid
title('LF8 DOWN Transition')
legend('L3MEC','L3LEC','LF8','LF2r')
yl = ylim();
line([0 0],yl,'color','k')
xlabel('Time lag (s)','fontsize',14)
ylabel('Amplitude (z)','fontsize',14)

