clear all
close all
%%
load F:\WC_Germany\overall_EC\overall_allcells_dir
addpath('F:\Code\WC_anal\general\')
addpath('F:\WC_Germany\Overall_EC\')
addpath('F:\Code\Chronux\spectral_analysis\continuous\')
addpath('F:\WC_Germany\hsmm_state_detection\\')
save_dir = 'F:\WC_Germany\overall_EC\allcells_lf8trig\';

drive_letter = 'F';
cd F:\WC_Germany\overall_EC

%%
dsf = 8;
Fsd = 2016/dsf;
niqf = 2016/2;
[b,a] = butter(2,[0.05/niqf 10/niqf]);

forwardlag = round(Fsd*1);
backwardlag = round(Fsd*1);
lags = -backwardlag:forwardlag;

for d = 1:length(sess_data)
    
    cdir = sess_data(d).directory;
    cdir(1) = 'F';
    disp(sprintf('session %d',d))
    cd(cdir);
    s_name = strcat(sess_data(d).region,'_l',sess_data(d).layer,'_',sess_data(d).name);
    load used_data lf2 lf3 lf4 lf5 lf8 
    
    lf8 = filtfilt(b,a,lf8);
    lf8 = downsample(lf8,dsf)/sess_data(d).gains(8);
    lf2 = filtfilt(b,a,lf2);
    lf2 = downsample(lf2,dsf)/sess_data(d).gains(2);
    lf3 = filtfilt(b,a,lf3);
    lf3 = downsample(lf3,dsf)/sess_data(d).gains(3);
    if ~exist('lf4','var') | sess_data(d).thom_elec==1
        lf4 = nan(size(lf3));
    else
        lf4 = filtfilt(b,a,lf4);
        lf4 = downsample(lf4,dsf)/sess_data(d).gains(4);
    end
    lf5 = filtfilt(b,a,lf5);
    lf5 = downsample(lf5,dsf)/sess_data(d).gains(5);
    
    lf2_r = lf2 - lf5;
    
    t_axis = (1:length(lf3))/Fsd;
    
    %% extract up and down transition times for MP and LF8
    load ec_hmm_state_seq8
    lf8_state_seq_c = hmm_bbstate_seq8;
    [new_lf8_seg_inds] = round(resample_uds_seg_inds(hmm8.UDS_segs,50.4,Fsd,length(t_axis)));
    
    lf8_state_seq = nan(size(t_axis));
    
    lf8_utrans = [];
    lf8_dtrans = [];
    for n = 1:hmm8.Nsegs
        lf8_state_seq(new_lf8_seg_inds(n,1):new_lf8_seg_inds(n,2)) = lf8_state_seq_c{n};
        cur_lf8_utrans = new_lf8_seg_inds(n,1) + find(lf8_state_seq_c{n}(1:end-1) == 1 & lf8_state_seq_c{n}(2:end) == 2);
        cur_lf8_dtrans = new_lf8_seg_inds(n,1) + find(lf8_state_seq_c{n}(1:end-1) == 2 & lf8_state_seq_c{n}(2:end) == 1);
        cur_lf8_dtrans(cur_lf8_dtrans < cur_lf8_utrans(1)) = [];
        cur_lf8_utrans(cur_lf8_utrans > cur_lf8_dtrans(end)) = [];
        lf8_utrans = [lf8_utrans; cur_lf8_utrans];
        lf8_dtrans = [lf8_dtrans; cur_lf8_dtrans];
    end
    

    %%  
    lf8_utrig_lf8(d,:) = get_event_trig_avg(lf8,lf8_utrans,forwardlag,backwardlag);
    lf8_dtrig_lf8(d,:) = get_event_trig_avg(lf8,lf8_dtrans,forwardlag,backwardlag);
    lf8_utrig_lf3(d,:) = get_event_trig_avg(lf3,lf8_utrans,forwardlag,backwardlag);
    lf8_dtrig_lf3(d,:) = get_event_trig_avg(lf3,lf8_dtrans,forwardlag,backwardlag);
    lf8_utrig_lf2(d,:) = get_event_trig_avg(lf2,lf8_utrans,forwardlag,backwardlag);
    lf8_dtrig_lf2(d,:) = get_event_trig_avg(lf2,lf8_dtrans,forwardlag,backwardlag);
    lf8_utrig_lf4(d,:) = get_event_trig_avg(lf4,lf8_utrans,forwardlag,backwardlag);
    lf8_dtrig_lf4(d,:) = get_event_trig_avg(lf4,lf8_dtrans,forwardlag,backwardlag);
     lf8_utrig_lf5(d,:) = get_event_trig_avg(lf5,lf8_utrans,forwardlag,backwardlag);
    lf8_dtrig_lf5(d,:) = get_event_trig_avg(lf5,lf8_dtrans,forwardlag,backwardlag);
    lf8_utrig_lf2_r(d,:) = get_event_trig_avg(lf2_r,lf8_utrans,forwardlag,backwardlag);
    lf8_dtrig_lf2_r(d,:) = get_event_trig_avg(lf2_r,lf8_dtrans,forwardlag,backwardlag);
 
    %%
        figure('visible','off')
    plot(lags/Fsd,lf8_utrig_lf8(d,:),'r'), hold on
    plot(lags/Fsd,lf8_utrig_lf2(d,:),'k')
    plot(lags/Fsd,lf8_utrig_lf2_r(d,:),'k--')
    plot(lags/Fsd,lf8_utrig_lf3(d,:),'g')
    plot(lags/Fsd,lf8_utrig_lf4(d,:),'c')
    plot(lags/Fsd,lf8_utrig_lf5(d,:),'b')
    legend('Lf8','LF2','LF2r','LF3','LF4','LF5')
    yl = ylim();
    line([0 0],yl,'color','k')
    savename = strcat(save_dir,'utrig_',s_name);
    print('-dpng',savename), close

        figure('visible','off')
    plot(lags/Fsd,lf8_dtrig_lf8(d,:),'r'), hold on
    plot(lags/Fsd,lf8_dtrig_lf2(d,:),'k')
    plot(lags/Fsd,lf8_dtrig_lf2_r(d,:),'k--')
    plot(lags/Fsd,lf8_dtrig_lf3(d,:),'g')
    plot(lags/Fsd,lf8_dtrig_lf4(d,:),'c')
    plot(lags/Fsd,lf8_dtrig_lf5(d,:),'b')
    legend('Lf8','LF2','LF2r','LF3','LF4','LF5')
    yl = ylim();
    line([0 0],yl,'color','k')
    savename = strcat(save_dir,'dtrig_',s_name);
    print('-dpng',savename), close
    
    clear lf2 lf2_r lf3 lf4 lf5 lf8
end

%%
cd F:\WC_Germany\overall_EC\
save overall_allcells_triglf3_data lags *_utrig_* *_dtrig_*

%%

% load overall_EC_coherence_data_gains
% mec = find_struct_field_vals(sess_data,'region','MEC');
% layer3 = find_struct_field_vals(sess_data,'layer','3');
% layer2 = find_struct_field_vals(sess_data,'layer','2');
% layer23 = find_struct_field_vals(sess_data,'layer','23');
% lec = find_struct_field_vals(sess_data,'region','LEC');
% l3mec = intersect(mec,layer3);
% l2mec = intersect(mec,layer2);
% l3lec = intersect(lec,layer3);
% l3mec(24:end) = [];
% l23mec = intersect(mec,layer23);
% l23mec = unique([l2mec l3mec l23mec]);
% 
% frange = find(f_i > 0.2 & f_i < 0.6);
% avg_P83 = mean(P83(:,frange),2);
% correct_lf3phase = find(avg_P83 > -0.4);
% l3mec_c = l3mec(avg_P83(l3mec) > -0.4);
% l3mec_w = l3mec(avg_P83(l3mec) < -0.4);
% l3lec_c = l3lec(avg_P83(l3lec) > -0.4);
% l3lec_w = l3lec(avg_P83(l3lec) < -0.4);

%%
% figure
% h = errorbar(lags/Fsd,nanmean(mp_utrig_mp(l3mec_c,:)),nanstd(mp_utrig_mp(l3mec_c,:))/sqrt(length(l3mec_c)));
% errorbar_tick(h,.01,'units');
% hold on
% h = errorbar(lags/Fsd,nanmean(lf8_utrig_lf8(l3mec_c,:)),nanstd(lf8_utrig_lf8(l3mec_c,:))/sqrt(length(l3mec_c)),'r');
% errorbar_tick(h,.01,'units');
% hold on
% h = errorbar(lags/Fsd,nanmean(mp_utrig_lf3(l3mec_c,:)),nanstd(mp_utrig_lf3(l3mec_c,:))/sqrt(length(l3mec_c)),'c');
% errorbar_tick(h,.01,'units');
% hold on
% h = errorbar(lags/Fsd,nanmean(lf8_utrig_lf3(l3mec_c,:)),nanstd(lf8_utrig_lf3(l3mec_c,:))/sqrt(length(l3mec_c)),'k');
% errorbar_tick(h,.01,'units');
% hold on
% xlim([-0.5 1])
% title('L3MEC Up Transition')
% legend('MP Up-trig MP','LF8 Up-trig Lf8','MP Up-trig LF3','LF8 Up-trig LF3')
% yl = ylim();
% line([0 0],yl,'color','k')
% xlabel('Time lag (s)')
% ylabel('Amplitude (z)')
% 
% 
% figure
% h = errorbar(lags/Fsd,nanmean(mp_dtrig_mp(l3mec_c,:)),nanstd(mp_dtrig_mp(l3mec_c,:))/sqrt(length(l3mec_c)));
% errorbar_tick(h,.01,'units');
% hold on
% h = errorbar(lags/Fsd,nanmean(lf8_dtrig_lf8(l3mec_c,:)),nanstd(lf8_dtrig_lf8(l3mec_c,:))/sqrt(length(l3mec_c)),'r');
% errorbar_tick(h,.01,'units');
% hold on
% h = errorbar(lags/Fsd,nanmean(mp_dtrig_lf3(l3mec_c,:)),nanstd(mp_dtrig_lf3(l3mec_c,:))/sqrt(length(l3mec_c)),'c');
% errorbar_tick(h,.01,'units');
% hold on
% h = errorbar(lags/Fsd,nanmean(lf8_dtrig_lf3(l3mec_c,:)),nanstd(lf8_dtrig_lf3(l3mec_c,:))/sqrt(length(l3mec_c)),'k');
% errorbar_tick(h,.01,'units');
% hold on
% xlim([-0.5 1])
% title('L3MEC Down Transition')
% legend('MP Down-trig MP','LF8 Down-trig Lf8','MP Down-trig LF3','LF8 Down-trig LF3')
% yl = ylim();
% line([0 0],yl,'color','k')
% xlabel('Time lag (s)')
% ylabel('Amplitude (z)')
% 
% 
% figure
% h = errorbar(lags/Fsd,nanmean(mp_utrig_mp(l3lec_c,:)),nanstd(mp_utrig_mp(l3lec_c,:))/sqrt(length(l3lec_c)));
% errorbar_tick(h,.01,'units');
% hold on
% h = errorbar(lags/Fsd,nanmean(lf8_utrig_lf8(l3lec_c,:)),nanstd(lf8_utrig_lf8(l3lec_c,:))/sqrt(length(l3lec_c)),'r');
% errorbar_tick(h,.01,'units');
% hold on
% h = errorbar(lags/Fsd,nanmean(mp_utrig_lf3(l3lec_c,:)),nanstd(mp_utrig_lf3(l3lec_c,:))/sqrt(length(l3lec_c)),'c');
% errorbar_tick(h,.01,'units');
% hold on
% h = errorbar(lags/Fsd,nanmean(lf8_utrig_lf3(l3lec_c,:)),nanstd(lf8_utrig_lf3(l3lec_c,:))/sqrt(length(l3lec_c)),'k');
% errorbar_tick(h,.01,'units');
% hold on
% xlim([-0.5 1])
% title('L3LEC Up Transition')
% legend('MP Up-trig MP','LF8 Up-trig Lf8','MP Up-trig LF3','LF8 Up-trig LF3')
% yl = ylim();
% line([0 0],yl,'color','k')
% xlabel('Time lag (s)')
% ylabel('Amplitude (z)')
% 
% 
% figure
% h = errorbar(lags/Fsd,nanmean(mp_dtrig_mp(l3lec_c,:)),nanstd(mp_dtrig_mp(l3lec_c,:))/sqrt(length(l3lec_c)));
% errorbar_tick(h,.01,'units');
% hold on
% h = errorbar(lags/Fsd,nanmean(lf8_dtrig_lf8(l3lec_c,:)),nanstd(lf8_dtrig_lf8(l3lec_c,:))/sqrt(length(l3lec_c)),'r');
% errorbar_tick(h,.01,'units');
% hold on
% h = errorbar(lags/Fsd,nanmean(mp_dtrig_lf3(l3lec_c,:)),nanstd(mp_dtrig_lf3(l3lec_c,:))/sqrt(length(l3lec_c)),'c');
% errorbar_tick(h,.01,'units');
% hold on
% h = errorbar(lags/Fsd,nanmean(lf8_dtrig_lf3(l3lec_c,:)),nanstd(lf8_dtrig_lf3(l3lec_c,:))/sqrt(length(l3lec_c)),'k');
% errorbar_tick(h,.01,'units');
% hold on
% xlim([-0.5 1])
% title('L3LEC Down Transition')
% legend('MP Down-trig MP','LF8 Down-trig Lf8','MP Down-trig LF3','LF8 Down-trig LF3')
% yl = ylim();
% line([0 0],yl,'color','k')
% xlabel('Time lag (s)')
% ylabel('Amplitude (z)')
% 
