clear all
close all

load C:\WC_Germany\overall_EC\overall_EC_dir
addpath('C:\WC_Germany\parietal_cortical_2010\')
addpath('C:\WC_Germany\hsmm_state_detection\')

used_data = [l3mec_p l3lec_p];
sess_data = sess_data(used_data);

drive_letter = 'C';

raw_Fs = 2016;
niqf = raw_Fs/2;
lcf = 0.05/niqf;
hcf = 40/niqf;
[b,a] = butter(2,[lcf hcf]);
dsf = 8;
Fsd = raw_Fs/dsf;

backlag = 5*Fsd;
forwardlag = 5*Fsd;
lags = (-backlag:forwardlag)/Fsd;

cd C:\WC_Germany\persistent_9_27_2010\
load pa_corresponding_lfp_revised_simp_new2

min_n_states = 1;
for d = 1:length(sess_data)
    
    cdir = sess_data(d).directory;
    cdir(1) = 'C';
    disp(sprintf('session %d',d))
    cd(cdir);
    s_name = strcat(sess_data(d).region,'_l',sess_data(d).layer,'_',sess_data(d).name);
    
    load ./used_data lf8 wcv_minus_spike lf3 lf5 lf2
        
    wcv_f = filtfilt(b,a,wcv_minus_spike);
    lf8_f = filtfilt(b,a,lf8);
    lf3_f = filtfilt(b,a,lf3);
%     lf2_f = filtfilt(b,a,lf2);
    
    down_w = downsample(wcv_f,dsf);
    down_8 = downsample(lf8_f,dsf);
    down_3 = downsample(lf3_f,dsf);
%     down_2 = downsample(lf2_f,dsf);
    
    down_w = zscore(down_w);
    down_8 = zscore(down_8);
    down_3 = zscore(down_3);
%     down_2 = zscore(down_2);
    
%     down_2hf = get_hf_features(lf2,raw_Fs,Fsd,[15 80],0.05);
%     down_3hf = get_hf_features(lf3,raw_Fs,Fsd,[15 80],0.05);
    down_3hf = get_hf_features(lf3,raw_Fs,Fsd,[15 80],0.05);
    down_2hf = get_hf_features(lf2,raw_Fs,Fsd,[15 80],0.05);
    
    
    load ./pa_hsmm_state_seq8_new2
    load ./pa_hsmm_state_seq_new2
    lfp_state_seq = hsmm_bbstate_seq8;
    mp_state_seq = hsmm_bbstate_seq;
    
    [new_seg_inds] = resample_uds_seg_inds(hmm8.UDS_segs,hmm8.Fs,252,length(down_8));
    [up_trans_inds8,down_trans_inds8] = compute_state_transitions_seg(new_seg_inds,lfp_state_seq);
    [up_trans_inds,down_trans_inds] = compute_state_transitions_seg(new_seg_inds,mp_state_seq);
    
    lf8_utrig_lf8(d,:) = get_event_trig_avg(down_8,up_trans_inds8,forwardlag,backlag);
    lf8_dtrig_lf8(d,:) = get_event_trig_avg(down_8,down_trans_inds8,forwardlag,backlag);
    lf8_utrig_mp(d,:) = get_event_trig_avg(down_w,up_trans_inds8,forwardlag,backlag);
    lf8_dtrig_mp(d,:) = get_event_trig_avg(down_w,down_trans_inds8,forwardlag,backlag);
    lf8_utrig_lf3(d,:) = get_event_trig_avg(down_3,up_trans_inds8,forwardlag,backlag);
    lf8_dtrig_lf3(d,:) = get_event_trig_avg(down_3,down_trans_inds8,forwardlag,backlag);
%     lf8_utrig_lf2(d,:) = get_event_trig_avg(down_2,up_trans_inds8,forwardlag,backlag);
%     lf8_dtrig_lf2(d,:) = get_event_trig_avg(down_2,down_trans_inds8,forwardlag,backlag);
    lf8_utrig_lf3h(d,:) = get_event_trig_avg(down_3hf,up_trans_inds8,forwardlag,backlag);
    lf8_dtrig_lf3h(d,:) = get_event_trig_avg(down_3hf,down_trans_inds8,forwardlag,backlag);
    lf8_utrig_lf2h(d,:) = get_event_trig_avg(down_2hf,up_trans_inds8,forwardlag,backlag);
    lf8_dtrig_lf2h(d,:) = get_event_trig_avg(down_2hf,down_trans_inds8,forwardlag,backlag);
%     lf8_utrig_lf2h(d,:) = get_event_trig_avg(down_2hf,up_trans_inds8,forwardlag,backlag);
%     lf8_dtrig_lf2h(d,:) = get_event_trig_avg(down_2hf,down_trans_inds8,forwardlag,backlag);

    mp_utrig_lf8(d,:) = get_event_trig_avg(down_8,up_trans_inds,forwardlag,backlag);
    mp_dtrig_lf8(d,:) = get_event_trig_avg(down_8,down_trans_inds,forwardlag,backlag);
    mp_utrig_mp(d,:) = get_event_trig_avg(down_w,up_trans_inds,forwardlag,backlag);
    mp_dtrig_mp(d,:) = get_event_trig_avg(down_w,down_trans_inds,forwardlag,backlag);
    mp_utrig_lf3(d,:) = get_event_trig_avg(down_3,up_trans_inds,forwardlag,backlag);
    mp_dtrig_lf3(d,:) = get_event_trig_avg(down_3,down_trans_inds,forwardlag,backlag);
%     mp_utrig_lf2(d,:) = get_event_trig_avg(down_2,up_trans_inds,forwardlag,backlag);
%     mp_dtrig_lf2(d,:) = get_event_trig_avg(down_2,down_trans_inds,forwardlag,backlag);
    mp_utrig_lf3h(d,:) = get_event_trig_avg(down_3hf,up_trans_inds,forwardlag,backlag);
    mp_dtrig_lf3h(d,:) = get_event_trig_avg(down_3hf,down_trans_inds,forwardlag,backlag);
    mp_utrig_lf2h(d,:) = get_event_trig_avg(down_2hf,up_trans_inds,forwardlag,backlag);
    mp_dtrig_lf2h(d,:) = get_event_trig_avg(down_2hf,down_trans_inds,forwardlag,backlag);
%     mp_utrig_lf2h(d,:) = get_event_trig_avg(down_2hf,up_trans_inds,forwardlag,backlag);
%     mp_dtrig_lf2h(d,:) = get_event_trig_avg(down_2hf,down_trans_inds,forwardlag,backlag);

    skipped_lfp_trans = [];
    for j = 1:length(up_trans_inds)
        skipped_lfp_trans = [skipped_lfp_trans mp_upskipped{d}.inds{j}];
    end
    non_skipped_lfp_trans = setdiff(1:length(down_trans_inds8),skipped_lfp_trans);
    n_uskip(d) = length(skipped_lfp_trans);

    dskipped_lfp_trans = [];
    for j = 1:length(up_trans_inds)
        dskipped_lfp_trans = [dskipped_lfp_trans mp_downskipped{d}.inds{j}];
    end
    non_dskipped_lfp_trans = setdiff(1:length(up_trans_inds8),dskipped_lfp_trans);
    n_dskip(d) = length(dskipped_lfp_trans);
    
    if length(dskipped_lfp_trans) >= min_n_states
        lf8_utrig_lf8_dsk(d,:) = get_event_trig_avg(down_8,up_trans_inds8(dskipped_lfp_trans),forwardlag,backlag);
        lf8_utrig_lf3_dsk(d,:) = get_event_trig_avg(down_3,up_trans_inds8(dskipped_lfp_trans),forwardlag,backlag);
%         lf8_utrig_lf2_dsk(d,:) = get_event_trig_avg(down_2,up_trans_inds8(dskipped_lfp_trans),forwardlag,backlag);
        lf8_utrig_lf3h_dsk(d,:) = get_event_trig_avg(down_3hf,up_trans_inds8(dskipped_lfp_trans),forwardlag,backlag);
        lf8_utrig_lf2h_dsk(d,:) = get_event_trig_avg(down_2hf,up_trans_inds8(dskipped_lfp_trans),forwardlag,backlag);
%         lf8_utrig_lf2h_dsk(d,:) = get_event_trig_avg(down_2hf,up_trans_inds8(dskipped_lfp_trans),forwardlag,backlag);
        lf8_utrig_mp_dsk(d,:) = get_event_trig_avg(down_w,up_trans_inds8(dskipped_lfp_trans),forwardlag,backlag);
    else
        lf8_utrig_lf8_dsk(d,:) = nan(1,forwardlag+backlag+1);
        lf8_utrig_lf3_dsk(d,:) = nan(1,forwardlag+backlag+1);
        lf8_utrig_mp_dsk(d,:) = nan(1,forwardlag+backlag+1);
%         lf8_utrig_lf2_dsk(d,:) = nan(1,forwardlag+backlag+1);
        lf8_utrig_lf3h_dsk(d,:) = nan(1,forwardlag+backlag+1);
        lf8_utrig_lf2h_dsk(d,:) = nan(1,forwardlag+backlag+1);
%         lf8_utrig_lf2h_dsk(d,:) = nan(1,forwardlag+backlag+1);
    end
    lf8_utrig_lf8_ndsk(d,:) = get_event_trig_avg(down_8,up_trans_inds8(non_dskipped_lfp_trans),forwardlag,backlag);
    lf8_utrig_mp_ndsk(d,:) = get_event_trig_avg(down_w,up_trans_inds8(non_dskipped_lfp_trans),forwardlag,backlag);
    lf8_utrig_lf3_ndsk(d,:) = get_event_trig_avg(down_3,up_trans_inds8(non_dskipped_lfp_trans),forwardlag,backlag);
%     lf8_utrig_lf2_ndsk(d,:) = get_event_trig_avg(down_2,up_trans_inds8(non_dskipped_lfp_trans),forwardlag,backlag);
    lf8_utrig_lf3h_ndsk(d,:) = get_event_trig_avg(down_3hf,up_trans_inds8(non_dskipped_lfp_trans),forwardlag,backlag);
    lf8_utrig_lf2h_ndsk(d,:) = get_event_trig_avg(down_2hf,up_trans_inds8(non_dskipped_lfp_trans),forwardlag,backlag);
%     lf8_utrig_lf2h_ndsk(d,:) = get_event_trig_avg(down_2hf,up_trans_inds8(non_dskipped_lfp_trans),forwardlag,backlag);

    if length(skipped_lfp_trans) >= min_n_states
        lf8_dtrig_lf8_sk(d,:) = get_event_trig_avg(down_8,down_trans_inds8(skipped_lfp_trans),forwardlag,backlag);
        lf8_dtrig_lf3_sk(d,:) = get_event_trig_avg(down_3,down_trans_inds8(skipped_lfp_trans),forwardlag,backlag);
%         lf8_dtrig_lf2_sk(d,:) = get_event_trig_avg(down_2,down_trans_inds8(skipped_lfp_trans),forwardlag,backlag);
        lf8_dtrig_mp_sk(d,:) = get_event_trig_avg(down_w,down_trans_inds8(skipped_lfp_trans),forwardlag,backlag);
        lf8_dtrig_lf3h_sk(d,:) = get_event_trig_avg(down_3hf,down_trans_inds8(skipped_lfp_trans),forwardlag,backlag);
        lf8_dtrig_lf2h_sk(d,:) = get_event_trig_avg(down_2hf,down_trans_inds8(skipped_lfp_trans),forwardlag,backlag);
%         lf8_dtrig_lf2h_sk(d,:) = get_event_trig_avg(down_2hf,down_trans_inds8(skipped_lfp_trans),forwardlag,backlag);
    else
        lf8_dtrig_lf8_sk(d,:) = nan;
        lf8_dtrig_mp_sk(d,:) = nan;
        lf8_dtrig_lf3_sk(d,:) = nan;
%         lf8_dtrig_lf2_sk(d,:) = nan;
        lf8_dtrig_lf3h_sk(d,:) = nan;
        lf8_dtrig_lf2h_sk(d,:) = nan;
%         lf8_dtrig_lf2h_sk(d,:) = nan;
    end
    lf8_dtrig_lf8_nsk(d,:) = get_event_trig_avg(down_8,down_trans_inds8(non_skipped_lfp_trans),forwardlag,backlag);
    lf8_dtrig_mp_nsk(d,:) = get_event_trig_avg(down_w,down_trans_inds8(non_skipped_lfp_trans),forwardlag,backlag);
    lf8_dtrig_lf3_nsk(d,:) = get_event_trig_avg(down_3,down_trans_inds8(non_skipped_lfp_trans),forwardlag,backlag);
%     lf8_dtrig_lf2_nsk(d,:) = get_event_trig_avg(down_2,down_trans_inds8(non_skipped_lfp_trans),forwardlag,backlag);
    lf8_dtrig_lf3h_nsk(d,:) = get_event_trig_avg(down_3hf,down_trans_inds8(non_skipped_lfp_trans),forwardlag,backlag);
    lf8_dtrig_lf2h_nsk(d,:) = get_event_trig_avg(down_2hf,down_trans_inds8(non_skipped_lfp_trans),forwardlag,backlag);
%     lf8_dtrig_lf2h_nsk(d,:) = get_event_trig_avg(down_2hf,down_trans_inds8(non_skipped_lfp_trans),forwardlag,backlag);
    
    skipped_lfp_trans(skipped_lfp_trans == length(up_trans_inds8)) = [];
    non_skipped_lfp_trans(non_skipped_lfp_trans == length(up_trans_inds8)) = [];
%     dskipped_lfp_trans(dskipped_lfp_trans == length(up_trans_inds8)) = [];
%     non_dskipped_lfp_trans(non_dskipped_lfp_trans == length(up_trans_inds8)) = [];
    if length(skipped_lfp_trans) >= min_n_states
        lf8_utrig_lf8_sk(d,:) = get_event_trig_avg(down_8,up_trans_inds8(skipped_lfp_trans+1),forwardlag,backlag);
        lf8_utrig_mp_sk(d,:) = get_event_trig_avg(down_w,up_trans_inds8(skipped_lfp_trans+1),forwardlag,backlag);
        lf8_utrig_lf3_sk(d,:) = get_event_trig_avg(down_3,up_trans_inds8(skipped_lfp_trans+1),forwardlag,backlag);
%         lf8_utrig_lf2_sk(d,:) = get_event_trig_avg(down_2,up_trans_inds8(skipped_lfp_trans+1),forwardlag,backlag);
        lf8_utrig_lf3h_sk(d,:) = get_event_trig_avg(down_3hf,up_trans_inds8(skipped_lfp_trans+1),forwardlag,backlag);
        lf8_utrig_lf2h_sk(d,:) = get_event_trig_avg(down_2hf,up_trans_inds8(skipped_lfp_trans+1),forwardlag,backlag);
%         lf8_utrig_lf2h_sk(d,:) = get_event_trig_avg(down_2hf,up_trans_inds8(skipped_lfp_trans+1),forwardlag,backlag);
    else
        lf8_utrig_lf8_sk(d,:) = nan;
        lf8_utrig_mp_sk(d,:) = nan;
        lf8_utrig_lf3_sk(d,:) = nan;
%         lf8_utrig_lf2_sk(d,:) = nan;
        lf8_utrig_lf3h_sk(d,:) = nan;
%         lf8_utrig_lf2h_sk(d,:) = nan;
    end
    lf8_utrig_lf8_nsk(d,:) = get_event_trig_avg(down_8,up_trans_inds8(non_skipped_lfp_trans+1),forwardlag,backlag);
    lf8_utrig_mp_nsk(d,:) = get_event_trig_avg(down_w,up_trans_inds8(non_skipped_lfp_trans+1),forwardlag,backlag);
    lf8_utrig_lf3_nsk(d,:) = get_event_trig_avg(down_3,up_trans_inds8(non_skipped_lfp_trans+1),forwardlag,backlag);
%     lf8_utrig_lf2_nsk(d,:) = get_event_trig_avg(down_2,up_trans_inds8(non_skipped_lfp_trans+1),forwardlag,backlag);
    lf8_utrig_lf3h_nsk(d,:) = get_event_trig_avg(down_3hf,up_trans_inds8(non_skipped_lfp_trans+1),forwardlag,backlag);
    lf8_utrig_lf2h_nsk(d,:) = get_event_trig_avg(down_2hf,up_trans_inds8(non_skipped_lfp_trans+1),forwardlag,backlag);
%     lf8_utrig_lf2h_nsk(d,:) = get_event_trig_avg(down_2hf,up_trans_inds8(non_skipped_lfp_trans+1),forwardlag,backlag);

    if length(dskipped_lfp_trans) >= min_n_states
        lf8_dtrig_lf8_dsk(d,:) = get_event_trig_avg(down_8,down_trans_inds8(dskipped_lfp_trans),forwardlag,backlag);
        lf8_dtrig_mp_dsk(d,:) = get_event_trig_avg(down_w,down_trans_inds8(dskipped_lfp_trans),forwardlag,backlag);
        lf8_dtrig_lf3_dsk(d,:) = get_event_trig_avg(down_3,down_trans_inds8(dskipped_lfp_trans),forwardlag,backlag);
%         lf8_dtrig_lf2_dsk(d,:) = get_event_trig_avg(down_2,down_trans_inds8(dskipped_lfp_trans),forwardlag,backlag);
        lf8_dtrig_lf3h_dsk(d,:) = get_event_trig_avg(down_3hf,down_trans_inds8(dskipped_lfp_trans),forwardlag,backlag);
        lf8_dtrig_lf2h_dsk(d,:) = get_event_trig_avg(down_2hf,down_trans_inds8(dskipped_lfp_trans),forwardlag,backlag);
%         lf8_dtrig_lf2h_dsk(d,:) = get_event_trig_avg(down_2hf,down_trans_inds8(dskipped_lfp_trans),forwardlag,backlag);
    else
        lf8_dtrig_lf8_dsk(d,:) = nan;
        lf8_dtrig_mp_dsk(d,:) = nan;
        lf8_dtrig_lf3_dsk(d,:) = nan;
%         lf8_dtrig_lf2_dsk(d,:) = nan;
        lf8_dtrig_lf3h_dsk(d,:) = nan;
        lf8_dtrig_lf2h_dsk(d,:) = nan;
%         lf8_dtrig_lf2h_dsk(d,:) = nan;
    end
    lf8_dtrig_lf8_ndsk(d,:) = get_event_trig_avg(down_8,down_trans_inds8(non_dskipped_lfp_trans),forwardlag,backlag);
    lf8_dtrig_mp_ndsk(d,:) = get_event_trig_avg(down_w,down_trans_inds8(non_dskipped_lfp_trans),forwardlag,backlag);
    lf8_dtrig_lf3_ndsk(d,:) = get_event_trig_avg(down_3,down_trans_inds8(non_dskipped_lfp_trans),forwardlag,backlag);
%     lf8_dtrig_lf2_ndsk(d,:) = get_event_trig_avg(down_2,down_trans_inds8(non_dskipped_lfp_trans),forwardlag,backlag);
    lf8_dtrig_lf3h_ndsk(d,:) = get_event_trig_avg(down_3hf,down_trans_inds8(non_dskipped_lfp_trans),forwardlag,backlag);
    lf8_dtrig_lf2h_ndsk(d,:) = get_event_trig_avg(down_2hf,down_trans_inds8(non_dskipped_lfp_trans),forwardlag,backlag);
%     lf8_dtrig_lf2h_ndsk(d,:) = get_event_trig_avg(down_2hf,down_trans_inds8(non_dskipped_lfp_trans),forwardlag,backlag);

end

clear *mat lf8_f lf8 lf3 lf3_f

cd C:\WC_Germany\persistent_9_27_2010\
save lf8_trig_avg_data_new3_lf8hftest2 lags lf8* mp* n_*

%%
mec = 1:22;
lec = 23:36;
bad_lf3 = [7 9 11 12 13 16 17 25 31 32 35];
good_lf3 = setdiff(1:36,bad_lf3);
good_lf3_mec = good_lf3(ismember(good_lf3,mec));
good_lf3_lec = good_lf3(ismember(good_lf3,lec));
%%

figure
h = errorbar(lags,nanmean(lf8_utrig_lf8),nanstd(lf8_utrig_lf8)/sqrt(36),'r')
hold on
errorbar_tick(h,.001,'units');
h = errorbar(lags,nanmean(lf8_utrig_lf3(good_lf3,:)),nanstd(lf8_utrig_lf3(good_lf3,:))/sqrt(length(good_lf3)),'k')
errorbar_tick(h,.001,'units');
h = errorbar(lags,nanmean(lf8_utrig_mp(mec,:)),nanstd(lf8_utrig_mp(mec,:))/sqrt(length(mec)),'b')
errorbar_tick(h,.001,'units');
h = errorbar(lags,nanmean(lf8_utrig_mp(lec,:)),nanstd(lf8_utrig_mp(lec,:))/sqrt(length(lec)),'g')
errorbar_tick(h,.001,'units');
% h = errorbar(lags,nanmean(lf8_utrig_lf2(good_lf3,:)),nanstd(lf8_utrig_lf2(good_lf3,:))/sqrt(length(good_lf3)),'c')
% errorbar_tick(h,.001,'units');
h = errorbar(lags,nanmean(lf8_utrig_lf3h(good_lf3,:)),nanstd(lf8_utrig_lf3h(good_lf3,:))/sqrt(length(good_lf3)),'y')
errorbar_tick(h,.001,'units');
h = errorbar(lags,nanmean(lf8_utrig_lf2h(good_lf3,:)),nanstd(lf8_utrig_lf2h(good_lf3,:))/sqrt(length(good_lf3)),'c')
errorbar_tick(h,.001,'units');

xlim([-1 2])
% legend('LF8-tr LF8','LF8-tr LF3','LF8-tr MEC','LF8-tr LEC')
yl = ylim;
line([0 0],yl,'color','k')

%%
figure
h = errorbar(lags,nanmean(lf8_dtrig_lf8),nanstd(lf8_dtrig_lf8)/sqrt(36),'r')
hold on
errorbar_tick(h,.001,'units');
h = errorbar(lags,nanmean(lf8_dtrig_mp(mec,:)),nanstd(lf8_dtrig_mp(mec,:))/sqrt(length(mec)),'b')
errorbar_tick(h,.001,'units');
h = errorbar(lags,nanmean(lf8_dtrig_mp(lec,:)),nanstd(lf8_dtrig_mp(lec,:))/sqrt(length(lec)),'g')
errorbar_tick(h,.001,'units');
h = errorbar(lags,nanmean(lf8_dtrig_lf3(good_lf3,:)),nanstd(lf8_dtrig_lf3(good_lf3,:))/sqrt(length(good_lf3)),'k')
errorbar_tick(h,.001,'units');
% h = errorbar(lags,nanmean(lf8_dtrig_lf2(good_lf3,:)),nanstd(lf8_dtrig_lf2(good_lf3,:))/sqrt(length(good_lf3)),'c')
% errorbar_tick(h,.001,'units');
h = errorbar(lags,nanmean(lf8_dtrig_lf3h(good_lf3,:)),nanstd(lf8_dtrig_lf3h(good_lf3,:))/sqrt(length(good_lf3)),'y')
errorbar_tick(h,.001,'units');
h = errorbar(lags,nanmean(lf8_dtrig_lf2h(good_lf3,:)),nanstd(lf8_dtrig_lf2h(good_lf3,:))/sqrt(length(good_lf3)),'c')
errorbar_tick(h,.001,'units');

xlim([-1 2])
% legend('LF8-tr LF8','LF8-tr LF3','LF8-tr MEC','LF8-tr LEC')
yl = ylim;
line([0 0],yl,'color','k')

%%
figure
h = errorbar(lags,nanmean(mp_utrig_lf8),nanstd(mp_utrig_lf8)/sqrt(36),'r')
errorbar_tick(h,.001,'units');
hold on
h = errorbar(lags,nanmean(mp_utrig_mp(mec,:)),nanstd(mp_utrig_mp(mec,:))/sqrt(length(mec)),'b')
errorbar_tick(h,.001,'units');
h = errorbar(lags,nanmean(mp_utrig_mp(lec,:)),nanstd(mp_utrig_mp(lec,:))/sqrt(length(lec)),'g')
errorbar_tick(h,.001,'units');
% errorbar_tick(h,.001,'units');
h = errorbar(lags,nanmean(mp_utrig_lf3(good_lf3_mec,:)),nanstd(mp_utrig_lf3(good_lf3_mec,:))/sqrt(length(good_lf3_mec)),'k')
errorbar_tick(h,.001,'units');
h = errorbar(lags,nanmean(mp_utrig_lf2h(good_lf3_lec,:)),nanstd(mp_utrig_lf2h(good_lf3_lec,:))/sqrt(length(good_lf3_lec)),'c')
errorbar_tick(h,.001,'units');
h = errorbar(lags,nanmean(mp_utrig_lf3h(good_lf3_lec,:)),nanstd(mp_utrig_lf3h(good_lf3_lec,:))/sqrt(length(good_lf3_lec)),'y')
errorbar_tick(h,.001,'units');
xlim([-2 3])
% legend('MEC-tr LF3','LEC-tr LF3','MEC-tr MP','LEC-tr MP')
yl = ylim;
line([0 0],yl,'color','k')

%%
figure
h = errorbar(lags,nanmean(mp_dtrig_lf8),nanstd(mp_dtrig_lf8)/sqrt(36),'r')
errorbar_tick(h,.001,'units');
hold on
h = errorbar(lags,nanmean(mp_dtrig_mp(mec,:)),nanstd(mp_dtrig_mp(mec,:))/sqrt(length(mec)),'b')
errorbar_tick(h,.001,'units');
h = errorbar(lags,nanmean(mp_dtrig_mp(lec,:)),nanstd(mp_dtrig_mp(lec,:))/sqrt(length(lec)),'g')
errorbar_tick(h,.001,'units');
h = errorbar(lags,nanmean(mp_dtrig_lf3(good_lf3_mec,:)),nanstd(mp_dtrig_lf3(good_lf3_mec,:))/sqrt(length(good_lf3_mec)),'k')
errorbar_tick(h,.001,'units');
h = errorbar(lags,nanmean(mp_dtrig_lf2(good_lf3_lec,:)),nanstd(mp_dtrig_lf3(good_lf3_lec,:))/sqrt(length(good_lf3_lec)),'c')
errorbar_tick(h,.001,'units');
h = errorbar(lags,nanmean(mp_dtrig_lf3h(good_lf3_mec,:)),nanstd(mp_dtrig_lf3h(good_lf3_mec,:))/sqrt(length(good_lf3_mec)),'y')
errorbar_tick(h,.001,'units');
% h = errorbar(lags,nanmean(mp_dtrig_lf2(good_lf3_mec,:)),nanstd(mp_dtrig_lf2(good_lf3_mec,:))/sqrt(length(good_lf3_mec)),'y')
% errorbar_tick(h,.001,'units');
xlim([-2 3])
% legend('MEC-tr LF3','LEC-tr LF3','MEC-tr MP','LEC-tr MP')
yl = ylim;
line([0 0],yl,'color','k')

figure
h = errorbar(lags,nanmean(lf8_utrig_lf3(good_lf3,:)),nanstd(lf8_utrig_lf3(good_lf3,:))/sqrt(length(good_lf3)),'r')
hold on
errorbar_tick(h,.001,'units');
h = errorbar(lags,nanmean(mp_utrig_lf3(good_lf3_mec,:)),nanstd(mp_utrig_lf3(good_lf3_mec,:))/sqrt(length(good_lf3_mec)),'b')
errorbar_tick(h,.001,'units');
h = errorbar(lags,nanmean(mp_utrig_lf3(good_lf3_lec,:)),nanstd(mp_utrig_lf3(good_lf3_lec,:))/sqrt(length(good_lf3_lec)),'g')
errorbar_tick(h,.001,'units');
h = errorbar(lags,nanmean(mp_utrig_lf3h(good_lf3_mec,:)),nanstd(mp_utrig_lf3h(good_lf3_mec,:))/sqrt(length(good_lf3_mec)),'y')
errorbar_tick(h,.001,'units');
% h = errorbar(lags,nanmean(mp_utrig_lf3h(good_lf3_lec,:)),nanstd(mp_utrig_lf3h(good_lf3_lec,:))/sqrt(length(good_lf3_lec)),'c')
% errorbar_tick(h,.001,'units');
xlim([-1 2])
legend('LF8-tr LF3','MEC-tr LF3','LEC-tr LF3')
yl = ylim;
line([0 0],yl,'color','k')

figure
h = errorbar(lags,nanmean(lf8_dtrig_lf3(good_lf3,:)),nanstd(lf8_dtrig_lf3(good_lf3,:))/sqrt(length(good_lf3)),'r')
hold on
errorbar_tick(h,.001,'units');
h = errorbar(lags,nanmean(mp_dtrig_lf3(good_lf3_mec,:)),nanstd(mp_dtrig_lf3(good_lf3_mec,:))/sqrt(length(good_lf3_mec)),'b')
errorbar_tick(h,.001,'units');
h = errorbar(lags,nanmean(mp_dtrig_lf3h(good_lf3_mec,:)),nanstd(mp_dtrig_lf3h(good_lf3_mec,:))/sqrt(length(good_lf3_mec)),'y')
errorbar_tick(h,.001,'units');
h = errorbar(lags,nanmean(mp_dtrig_lf3(good_lf3_lec,:)),nanstd(mp_dtrig_lf3(good_lf3_lec,:))/sqrt(length(good_lf3_lec)),'g')
errorbar_tick(h,.001,'units');
xlim([-1 2])
legend('LF8-tr LF3','MEC-tr LF3','LEC-tr LF3')
yl = ylim;
line([0 0],yl,'color','k')

%%
min_states = 10;
gu_mec = good_lf3_mec(n_uskip(good_lf3_mec) >= min_states);
gu_lec = good_lf3_lec(n_uskip(good_lf3_lec) >= min_states);
gd_mec = good_lf3_mec(n_uskip(good_lf3_mec) >= min_states);
gd_lec = good_lf3_lec(n_uskip(good_lf3_lec) >= min_states);

figure
set(gca,'fontsize',12,'fontname','arial')
h = errorbar(lags,nanmean(lf8_dtrig_lf8),nanstd(lf8_dtrig_lf8)/sqrt(36),'r')
hold on
errorbar_tick(h,.001,'units');
h = errorbar(lags,nanmean(lf8_dtrig_lf8_sk),nanstd(lf8_dtrig_lf8_sk)/sqrt(36),'k')
errorbar_tick(h,.001,'units');
h = errorbar(lags,nanmean(lf8_dtrig_lf2h_sk(gu_mec,:)),nanstd(lf8_dtrig_lf2h_sk(gu_mec,:))/sqrt(length(gu_mec)),'g')
errorbar_tick(h,.001,'units');
h = errorbar(lags,nanmean(lf8_dtrig_lf2h_nsk(gu_mec,:)),nanstd(lf8_dtrig_lf2h_nsk(gu_mec,:))/sqrt(length(gu_mec)),'color',[0.2 0.6 0.2])
errorbar_tick(h,.001,'units');
h = errorbar(lags,nanmean(lf8_dtrig_mp_sk(gu_mec,:)),nanstd(lf8_dtrig_mp_sk(gu_mec,:))/sqrt(length(gu_mec)),'b')
errorbar_tick(h,.001,'units');
h = errorbar(lags,nanmean(lf8_dtrig_mp_nsk(gu_mec,:)),nanstd(lf8_dtrig_mp_nsk(gu_mec,:))/sqrt(length(gu_mec)),'c')
errorbar_tick(h,.001,'units');
xlim([-2 4])
legend('LF8-tr LF8','LF8-tr LF3 Pers','Lf8-tr LF3 NPers','LF8-tr MP Pers','LF8-tr MP Npers')
yl = ylim;
line([0 0],yl,'color','k')
xlabel('Time (s)','fontname','arial')
ylabel('Amplitude (z)','fontname','arial')

%%
figure
set(gca,'fontsize',12,'fontname','arial')
h = errorbar(lags,nanmean(lf8_dtrig_lf8),nanstd(lf8_dtrig_lf8)/sqrt(36),'r')
hold on
errorbar_tick(h,.001,'units');
h = errorbar(lags,nanmean(lf8_dtrig_lf3_sk(gu_lec,:)),nanstd(lf8_dtrig_lf8_sk(gu_lec,:))/sqrt(length(gu_lec)),'k')
errorbar_tick(h,.001,'units');
h = errorbar(lags,nanmean(lf8_dtrig_lf3_nsk(gu_lec,:)),nanstd(lf8_dtrig_lf3_nsk(gu_lec,:))/sqrt(length(gu_lec)),'color',[0.2 0.6 0.2])
errorbar_tick(h,.001,'units');
h = errorbar(lags,nanmean(lf8_dtrig_mp_sk(gu_lec,:)),nanstd(lf8_dtrig_mp_sk(gu_lec,:))/sqrt(length(gu_lec)),'b')
errorbar_tick(h,.001,'units');
h = errorbar(lags,nanmean(lf8_dtrig_mp_nsk(gu_lec,:)),nanstd(lf8_dtrig_mp_nsk(gu_lec,:))/sqrt(length(gu_lec)),'c')
errorbar_tick(h,.001,'units');
xlim([-2 4])
legend('LF8-tr LF8','LF8-tr LF3 Pers','Lf8-tr LF3 NPers','LF8-tr MP Pers','LF8-tr MP Npers')
yl = ylim;
line([0 0],yl,'color','k')
xlabel('Time (s)','fontname','arial')
ylabel('Amplitude (z)','fontname','arial')

figure
set(gca,'fontsize',12,'fontname','arial')
h = errorbar(lags,nanmean(lf8_utrig_lf8),nanstd(lf8_utrig_lf8)/sqrt(36),'r')
hold on
errorbar_tick(h,.001,'units');
h = errorbar(lags,nanmean(lf8_utrig_lf3h_sk(gu_mec,:)),nanstd(lf8_utrig_lf3h_sk(gu_mec,:))/sqrt(length(gu_mec)),'g')
errorbar_tick(h,.001,'units');
h = errorbar(lags,nanmean(lf8_utrig_lf3h_nsk(gu_mec,:)),nanstd(lf8_utrig_lf3h_nsk(gu_mec,:))/sqrt(length(gu_mec)),'color',[0.2 0.6 0.2])
errorbar_tick(h,.001,'units');
h = errorbar(lags,nanmean(lf8_utrig_mp_sk(gu_mec,:)),nanstd(lf8_utrig_mp_sk(gu_mec,:))/sqrt(length(gu_mec)),'b')
errorbar_tick(h,.001,'units');
h = errorbar(lags,nanmean(lf8_utrig_mp_nsk(gu_mec,:)),nanstd(lf8_utrig_mp_nsk(gu_mec,:))/sqrt(length(gu_mec)),'c')
errorbar_tick(h,.001,'units');
xlim([-3 4])
legend('LF8-tr LF8','LF8-tr LF3 Pers','Lf8-tr LF3 NPers','LF8-tr MP Pers','LF8-tr MP Npers')
yl = ylim;
line([0 0],yl,'color','k')
xlabel('Time (s)','fontname','arial')
ylabel('Amplitude (z)','fontname','arial')

figure
h = errorbar(lags,nanmean(lf8_utrig_lf8),nanstd(lf8_utrig_lf8)/sqrt(36),'r')
hold on
errorbar_tick(h,.001,'units');
h = errorbar(lags,nanmean(lf8_utrig_lf3_sk(gu_lec,:)),nanstd(lf8_utrig_lf8_sk(gu_lec,:))/sqrt(length(gu_lec)),'k')
errorbar_tick(h,.001,'units');
h = errorbar(lags,nanmean(lf8_utrig_lf3_nsk(gu_lec,:)),nanstd(lf8_utrig_lf3_nsk(gu_lec,:))/sqrt(length(gu_lec)),'color',[0.2 0.6 0.2])
errorbar_tick(h,.001,'units');
h = errorbar(lags,nanmean(lf8_utrig_mp_sk(gu_lec,:)),nanstd(lf8_utrig_mp_sk(gu_lec,:))/sqrt(length(gu_lec)),'b')
errorbar_tick(h,.001,'units');
h = errorbar(lags,nanmean(lf8_utrig_mp_nsk(gu_lec,:)),nanstd(lf8_utrig_mp_nsk(gu_lec,:))/sqrt(length(gu_lec)),'c')
errorbar_tick(h,.001,'units');
xlim([-3 4])
legend('LF8-tr LF8','LF8-tr LF3 Pers','Lf8-tr LF3 NPers','LF8-tr MP Pers','LF8-tr MP Npers')
yl = ylim;
line([0 0],yl,'color','k')

figure
set(gca,'fontsize',12,'fontname','arial')
h = errorbar(lags,nanmean(lf8_utrig_lf8),nanstd(lf8_utrig_lf8)/sqrt(36),'r')
hold on
errorbar_tick(h,.001,'units');
h = errorbar(lags,nanmean(lf8_utrig_lf3h_dsk(gd_mec,:)),nanstd(lf8_utrig_lf3h_dsk(gd_mec,:))/sqrt(length(gd_mec)),'g')
errorbar_tick(h,.001,'units');
h = errorbar(lags,nanmean(lf8_utrig_lf3h_ndsk(gd_mec,:)),nanstd(lf8_utrig_lf3h_ndsk(gd_mec,:))/sqrt(length(gd_mec)),'color',[0.2 0.6 0.2])
errorbar_tick(h,.001,'units');
h = errorbar(lags,nanmean(lf8_utrig_mp_dsk(gd_mec,:)),nanstd(lf8_utrig_mp_dsk(gd_mec,:))/sqrt(length(gd_mec)),'b')
errorbar_tick(h,.001,'units');
h = errorbar(lags,nanmean(lf8_utrig_mp_ndsk(gd_mec,:)),nanstd(lf8_utrig_mp_ndsk(gd_mec,:))/sqrt(length(gd_mec)),'c')
errorbar_tick(h,.001,'units');
xlim([-2 4])
legend('LF8-tr LF8','LF8-tr LF3 Pers','Lf8-tr LF3 NPers','LF8-tr MP Pers','LF8-tr MP Npers')
yl = ylim;
line([0 0],yl,'color','k')
xlabel('Time (s)','fontname','arial')
ylabel('Amplitude (z)','fontname','arial')

figure
h = errorbar(lags,nanmean(lf8_utrig_lf8),nanstd(lf8_utrig_lf8)/sqrt(36),'r')
hold on
errorbar_tick(h,.001,'units');
h = errorbar(lags,nanmean(lf8_utrig_lf3_dsk(gd_lec,:)),nanstd(lf8_utrig_lf8_dsk(gd_lec,:))/sqrt(length(gd_lec)),'k')
errorbar_tick(h,.001,'units');
h = errorbar(lags,nanmean(lf8_utrig_lf3_ndsk(gd_lec,:)),nanstd(lf8_utrig_lf3_ndsk(gd_lec,:))/sqrt(length(gd_lec)),'color',[0.2 0.6 0.2])
errorbar_tick(h,.001,'units');
h = errorbar(lags,nanmean(lf8_utrig_mp_dsk(gd_lec,:)),nanstd(lf8_utrig_mp_dsk(gd_lec,:))/sqrt(length(gd_lec)),'b')
errorbar_tick(h,.001,'units');
h = errorbar(lags,nanmean(lf8_utrig_mp_ndsk(gd_lec,:)),nanstd(lf8_utrig_mp_ndsk(gd_lec,:))/sqrt(length(gd_lec)),'c')
errorbar_tick(h,.001,'units');
xlim([-3 4])
legend('LF8-tr LF8','LF8-tr LF3 Pers','Lf8-tr LF3 NPers','LF8-tr MP Pers','LF8-tr MP Npers')
yl = ylim;
line([0 0],yl,'color','k')

figure
set(gca,'fontsize',12,'fontname','arial')
h = errorbar(lags,nanmean(lf8_dtrig_lf8),nanstd(lf8_dtrig_lf8)/sqrt(36),'r')
hold on
errorbar_tick(h,.001,'units');
h = errorbar(lags,nanmean(lf8_dtrig_lf3_dsk(gd_mec,:)),nanstd(lf8_dtrig_lf8_dsk(gd_mec,:))/sqrt(length(gd_mec)),'g')
errorbar_tick(h,.001,'units');
h = errorbar(lags,nanmean(lf8_dtrig_lf3_ndsk(gd_mec,:)),nanstd(lf8_dtrig_lf3_ndsk(gd_mec,:))/sqrt(length(gd_mec)),'color',[0.2 0.6 0.2])
errorbar_tick(h,.001,'units');
h = errorbar(lags,nanmean(lf8_dtrig_mp_dsk(gd_mec,:)),nanstd(lf8_dtrig_mp_dsk(gd_mec,:))/sqrt(length(gd_mec)),'b')
errorbar_tick(h,.001,'units');
h = errorbar(lags,nanmean(lf8_dtrig_mp_ndsk(gd_mec,:)),nanstd(lf8_dtrig_mp_ndsk(gd_mec,:))/sqrt(length(gd_mec)),'c')
errorbar_tick(h,.001,'units');
xlim([-3 4])
legend('LF8-tr LF8','LF8-tr LF3 Pers','Lf8-tr LF3 NPers','LF8-tr MP Pers','LF8-tr MP Npers')
yl = ylim;
line([0 0],yl,'color','k')
xlabel('Time (s)','fontname','arial')
ylabel('Amplitude (z)','fontname','arial')

figure
h = errorbar(lags,nanmean(lf8_utrig_lf8),nanstd(lf8_utrig_lf8)/sqrt(36),'r')
hold on
errorbar_tick(h,.001,'units');
h = errorbar(lags,nanmean(lf8_utrig_lf3_dsk(gd_lec,:)),nanstd(lf8_utrig_lf8_dsk(gd_lec,:))/sqrt(length(gd_lec)),'g')
errorbar_tick(h,.001,'units');
h = errorbar(lags,nanmean(lf8_utrig_lf3_ndsk(gd_lec,:)),nanstd(lf8_utrig_lf3_ndsk(gd_lec,:))/sqrt(length(gd_lec)),'color',[0.2 0.6 0.2])
errorbar_tick(h,.001,'units');
h = errorbar(lags,nanmean(lf8_utrig_mp_dsk(gd_lec,:)),nanstd(lf8_utrig_mp_dsk(gd_lec,:))/sqrt(length(gd_lec)),'b')
errorbar_tick(h,.001,'units');
h = errorbar(lags,nanmean(lf8_utrig_mp_ndsk(gd_lec,:)),nanstd(lf8_utrig_mp_ndsk(gd_lec,:))/sqrt(length(gd_lec)),'c')
errorbar_tick(h,.001,'units');
xlim([-3 4])
legend('LF8-tr LF8','LF8-tr LF3 Pers','Lf8-tr LF3 NPers','LF8-tr MP Pers','LF8-tr MP Npers')
yl = ylim;
line([0 0],yl,'color','k')

%%
d =14;
figure
plot(lags,lf8_dtrig_lf8(d,:),'r')
hold on
plot(lags,lf8_dtrig_lf2_sk(d,:),'k')
plot(lags,lf8_dtrig_lf2_nsk(d,:),'color',[0.2 0.6 0.2])
plot(lags,lf8_dtrig_mp_sk(d,:),'b')
plot(lags,lf8_dtrig_mp_nsk(d,:),'c')
xlim([-2 4])
legend('LF8-tr LF8','LF8-tr LF3 Pers','Lf8-tr LF3 NPers','LF8-tr MP Pers','LF8-tr MP Npers')
yl = ylim;
line([0 0],yl,'color','k')

figure
plot(lags,lf8_dtrig_lf3(d,:),'r')
hold on
plot(lags,mp_dtrig_lf3(d,:),'b')
xlim([-2 3])
yl = ylim;
line([0 0],yl,'color','k')
legend('LF8-tr LF3','MP-tr LF3')

figure
plot(lags,lf8_utrig_lf3(d,:),'r')
hold on
plot(lags,mp_utrig_lf3(d,:),'b')
xlim([-2 3])
yl = ylim;
line([0 0],yl,'color','k')
legend('LF8-tr LF3','MP-tr LF3')


%%
figure
h = errorbar(lags,nanmean(lf8_utrig_lf8),nanstd(lf8_utrig_lf8)/sqrt(36),'r')
hold on
errorbar_tick(h,.001,'units');
% h = errorbar(lags,nanmean(lf8_utrig_lf3(good_lf3,:)),nanstd(lf8_utrig_lf3(good_lf3,:))/sqrt(length(good_lf3)),'k')
% errorbar_tick(h,.001,'units');
h = errorbar(lags,nanmean(lf8_utrig_lf3(good_lf3,:)),nanstd(lf8_utrig_lf3(good_lf3,:))/sqrt(length(good_lf3)),'color',[0.7 0.2 0.1])
errorbar_tick(h,.001,'units');
h = errorbar(lags,nanmean(lf8_utrig_lf2(good_lf3,:)),nanstd(lf8_utrig_lf2(good_lf3,:))/sqrt(length(good_lf3)),'color',[0.2 0.7 0.1])
errorbar_tick(h,.001,'units');
h = errorbar(lags,nanmean(lf8_utrig_mp(mec,:)),nanstd(lf8_utrig_mp(mec,:))/sqrt(length(mec)),'b')
errorbar_tick(h,.001,'units');
% h = errorbar(lags,nanmean(lf8_utrig_mp(lec,:)),nanstd(lf8_utrig_mp(lec,:))/sqrt(length(lec)),'g')
% errorbar_tick(h,.001,'units');
% h = errorbar(lags,nanmean(lf8_utrig_lf2(good_lf3,:)),nanstd(lf8_utrig_lf2(good_lf3,:))/sqrt(length(good_lf3)),'c')
% errorbar_tick(h,.001,'units');
h = errorbar(lags,nanmean(lf8_utrig_lf3h(good_lf3,:)),nanstd(lf8_utrig_lf3h(good_lf3,:))/sqrt(length(good_lf3)),'g')
errorbar_tick(h,.001,'units');
xlim([-2 3])

figure
h = errorbar(lags,nanmean(lf8_dtrig_lf8),nanstd(lf8_dtrig_lf8)/sqrt(36),'r')
hold on
errorbar_tick(h,.001,'units');
% h = errorbar(lags,nanmean(lf8_dtrig_lf3(good_lf3,:)),nanstd(lf8_dtrig_lf3(good_lf3,:))/sqrt(length(good_lf3)),'k')
% errorbar_tick(h,.001,'units');
h = errorbar(lags,nanmean(lf8_dtrig_lf3(good_lf3,:)),nanstd(lf8_dtrig_lf3(good_lf3,:))/sqrt(length(good_lf3)),'color',[0.7 0.2 0.1])
errorbar_tick(h,.001,'units');
h = errorbar(lags,nanmean(lf8_dtrig_lf2(good_lf3,:)),nanstd(lf8_dtrig_lf2(good_lf3,:))/sqrt(length(good_lf3)),'color',[0.2 0.7 0.1])
errorbar_tick(h,.001,'units');
h = errorbar(lags,nanmean(lf8_dtrig_mp(mec,:)),nanstd(lf8_dtrig_mp(mec,:))/sqrt(length(mec)),'b')
errorbar_tick(h,.001,'units');
% h = errorbar(lags,nanmean(lf8_utrig_mp(lec,:)),nanstd(lf8_utrig_mp(lec,:))/sqrt(length(lec)),'g')
% errorbar_tick(h,.001,'units');
% h = errorbar(lags,nanmean(lf8_utrig_lf2(good_lf3,:)),nanstd(lf8_utrig_lf2(good_lf3,:))/sqrt(length(good_lf3)),'c')
% errorbar_tick(h,.001,'units');
h = errorbar(lags,nanmean(lf8_dtrig_lf3h(good_lf3,:)),nanstd(lf8_dtrig_lf3h(good_lf3,:))/sqrt(length(good_lf3)),'g')
errorbar_tick(h,.001,'units');
xlim([-2 3])

figure
h = errorbar(lags,nanmean(mp_utrig_lf8(mec,:)),nanstd(mp_utrig_lf8(mec,:))/sqrt(length(mec)),'r')
errorbar_tick(h,.001,'units');
hold on
h = errorbar(lags,nanmean(mp_utrig_mp(mec,:)),nanstd(mp_utrig_mp(mec,:))/sqrt(length(mec)),'b')
errorbar_tick(h,.001,'units');
% h = errorbar(lags,nanmean(mp_utrig_lf3(good_lf3_mec,:)),nanstd(mp_utrig_lf3(good_lf3_mec,:))/sqrt(length(good_lf3_mec)),'k')
% errorbar_tick(h,.001,'units');
h = errorbar(lags,nanmean(mp_utrig_lf3(good_lf3_mec,:)),nanstd(mp_utrig_lf3(good_lf3_mec,:))/sqrt(length(good_lf3_mec)),'color',[0.7 0.2 0.1])
errorbar_tick(h,.001,'units');
h = errorbar(lags,nanmean(mp_utrig_lf2(good_lf3_mec,:)),nanstd(mp_utrig_lf2(good_lf3_mec,:))/sqrt(length(good_lf3_mec)),'color',[0.2 0.7 0.1])
errorbar_tick(h,.001,'units');
h = errorbar(lags,nanmean(mp_utrig_lf3h(good_lf3_mec,:)),nanstd(mp_utrig_lf3h(good_lf3_mec,:))/sqrt(length(good_lf3_mec)),'c')
errorbar_tick(h,.001,'units');
xlim([-2 3])

figure
h = errorbar(lags,nanmean(mp_dtrig_lf8(mec,:)),nanstd(mp_dtrig_lf8(mec,:))/sqrt(length(mec)),'r')
errorbar_tick(h,.001,'units');
hold on
h = errorbar(lags,nanmean(mp_dtrig_mp(mec,:)),nanstd(mp_dtrig_mp(mec,:))/sqrt(length(mec)),'b')
errorbar_tick(h,.001,'units');
h = errorbar(lags,nanmean(mp_dtrig_lf3(good_lf3_mec,:)),nanstd(mp_dtrig_lf3(good_lf3_mec,:))/sqrt(length(good_lf3_mec)),'color',[0.7 0.2 0.1])
errorbar_tick(h,.001,'units');
h = errorbar(lags,nanmean(mp_dtrig_lf2(good_lf3_mec,:)),nanstd(mp_dtrig_lf2(good_lf3_mec,:))/sqrt(length(good_lf3_mec)),'color',[0.2 0.7 0.1])
errorbar_tick(h,.001,'units');
h = errorbar(lags,nanmean(mp_dtrig_lf3h(good_lf3_mec,:)),nanstd(mp_dtrig_lf3h(good_lf3_mec,:))/sqrt(length(good_lf3_mec)),'g')
errorbar_tick(h,.001,'units');
xlim([-2 3])