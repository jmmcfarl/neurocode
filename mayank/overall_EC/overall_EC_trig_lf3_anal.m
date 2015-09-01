clear all
close all
%%
load E:\WC_Germany\overall_EC\overall_allcells_dir.mat
addpath('E:\Code\WC_anal\general\')
addpath('E:\WC_Germany\Overall_EC\')
addpath('E:\Code\Chronux\spectral_analysis\continuous\')
addpath('E:\WC_Germany\hsmm_state_detection\\')

drive_letter = 'E';
cd E:\WC_Germany\overall_EC
load corresponding_lfp_state_data
%%
raw_Fs = 2016;
dsf = 16;
Fsd = raw_Fs/dsf;
niqf = raw_Fs/2;

dsf_fac = round(252/Fsd);
amp_grid = linspace(-4,4,200);
forwardlag = round(Fsd*1);
backwardlag = round(Fsd*1);
lags = -backwardlag:forwardlag;

layer3 = find_struct_field_vals(sess_data,'layer','3');

for d = layer3
    
    cdir = sess_data(d).directory;
    cdir(1) = 'E';
    disp(sprintf('session %d',d))
    cd(cdir);
    s_name = strcat(sess_data(d).region,'_l',sess_data(d).layer,'_',sess_data(d).name);
    load used_data lf2 lf5 lf8 wcv_minus_spike
    load spike_time_jmm
    spk_id = round(spkid/dsf);
    
    lf2_r = lf2/sess_data(d).gains(2) - lf5/sess_data(d).gains(5);
    
    wcv = get_lf_features(wcv_minus_spike,raw_Fs,Fsd,[0.05 40]);
    lf8 = get_lf_features(lf8,raw_Fs,Fsd,[0.05 40]);
    lf2_rlf = get_lf_features(lf2_r,raw_Fs,Fsd,[0.05 40]);    
    lf2_hf = get_hf_features(lf2,raw_Fs,Fsd,[15 80],0.05);
    lf2_rhf = get_hf_features(lf2_r,raw_Fs,Fsd,[15 80],0.05);
    
    t_axis = (1:length(lf8))/Fsd;
    
    bin_spk = histc(spk_id,1:length(lf8));
    
    
    %% extract up and down transition times for MP and LF8
    load ec_hmm_state_seq
    load ec_hmm_state_seq8
    mp_state_seq_c =  hmm_bbstate_seq;
    lf8_state_seq_c = hmm_bbstate_seq8;
    [new_mp_seg_inds] = round(resample_uds_seg_inds(hmm.UDS_segs,50.4,Fsd,length(t_axis)));
    [new_lf8_seg_inds] = round(resample_uds_seg_inds(hmm8.UDS_segs,50.4,Fsd,length(t_axis)));
    
    for n = 1:hmm.Nsegs
        mp_state_seq_c{n} = downsample(hmm_bbstate_seq{n},dsf_fac);
    end
    for n = 1:hmm8.Nsegs
        lf8_state_seq_c{n} = downsample(hmm_bbstate_seq8{n},dsf_fac);
    end
    
    mp_state_seq = nan(size(t_axis));
    lf8_state_seq = nan(size(t_axis));
    mp_utrans = [];
    mp_dtrans = [];
    for n = 1:hmm.Nsegs
        mp_state_seq(new_mp_seg_inds(n,1):new_mp_seg_inds(n,1)+length(mp_state_seq_c{n})-1) = mp_state_seq_c{n};
        cur_mp_utrans = new_mp_seg_inds(n,1) + find(mp_state_seq_c{n}(1:end-1) == 1 & mp_state_seq_c{n}(2:end) == 2);
        cur_mp_dtrans = new_mp_seg_inds(n,1) + find(mp_state_seq_c{n}(1:end-1) == 2 & mp_state_seq_c{n}(2:end) == 1);
        cur_mp_dtrans(cur_mp_dtrans < cur_mp_utrans(1)) = [];
        cur_mp_utrans(cur_mp_utrans > cur_mp_dtrans(end)) = [];
        mp_utrans = [mp_utrans; cur_mp_utrans];
        mp_dtrans = [mp_dtrans; cur_mp_dtrans];
    end
    
    lf8_utrans = [];
    lf8_dtrans = [];
    for n = 1:hmm8.Nsegs
        lf8_state_seq(new_lf8_seg_inds(n,1):new_lf8_seg_inds(n,1)+length(lf8_state_seq_c{n})-1) = lf8_state_seq_c{n};
        cur_lf8_utrans = new_lf8_seg_inds(n,1) + find(lf8_state_seq_c{n}(1:end-1) == 1 & lf8_state_seq_c{n}(2:end) == 2);
        cur_lf8_dtrans = new_lf8_seg_inds(n,1) + find(lf8_state_seq_c{n}(1:end-1) == 2 & lf8_state_seq_c{n}(2:end) == 1);
        cur_lf8_dtrans(cur_lf8_dtrans < cur_lf8_utrans(1)) = [];
        cur_lf8_utrans(cur_lf8_utrans > cur_lf8_dtrans(end)) = [];
        lf8_utrans = [lf8_utrans; cur_lf8_utrans];
        lf8_dtrans = [lf8_dtrans; cur_lf8_dtrans];
    end
    

    %%
%     used_data = find(~isnan(lf8_state_seq) & ~isnan(mp_state_seq));    
%     used_data = find(~isnan(lf8_state_seq));    
%     p_lf8up(d) = sum(lf8_state_seq(used_data)==2)/length(used_data);
%     p_lf8down(d) = sum(lf8_state_seq(used_data)==1)/length(used_data);
%     p_mpup(d) = sum(mp_state_seq(used_data)==2)/length(used_data);
%     p_mpdown(d) = sum(mp_state_seq(used_data)==1)/length(used_data);
%     p_lf8up_cond_mpdown(d) = sum(lf8_state_seq(used_data)==2 & mp_state_seq(used_data)==1)/length(used_data)/p_mpdown(d);
%     p_lf8down_cond_mpdown(d) = sum(lf8_state_seq(used_data)==1 & mp_state_seq(used_data)==1)/length(used_data)/p_mpdown(d);
%     p_lf8up_cond_mpup(d) = sum(lf8_state_seq(used_data)==2 & mp_state_seq(used_data)==2)/length(used_data)/p_mpup(d);
%     p_lf8down_cond_mpup(d) = sum(lf8_state_seq(used_data)==1 & mp_state_seq(used_data)==2)/length(used_data)/p_mpup(d);
%     p_mpup_cond_lf8down(d) = sum(mp_state_seq(used_data)==2 & lf8_state_seq(used_data)==1)/length(used_data)/p_lf8down(d);
%     p_mpdown_cond_lf8down(d) = sum(mp_state_seq(used_data)==1 & lf8_state_seq(used_data)==1)/length(used_data)/p_lf8down(d);
%     p_mpup_cond_lf8up(d) = sum(mp_state_seq(used_data)==2 & lf8_state_seq(used_data)==2)/length(used_data)/p_lf8up(d);
%     p_mpdown_cond_lf8up(d) = sum(mp_state_seq(used_data)==1 & lf8_state_seq(used_data)==2)/length(used_data)/p_lf8up(d);    
    
%     lf8_utrans_cond_mup = lf8_utrans;
%     lf8_utrans_cond_mup(mp_state_seq(lf8_utrans_cond_mup) == 1) = [];
%     lf8_utrans_cond_mdown = lf8_utrans;
%     lf8_utrans_cond_mdown(mp_state_seq(lf8_utrans_cond_mdown) == 2) = [];
%     lf8_dtrans_cond_mup = lf8_dtrans;
%     lf8_dtrans_cond_mup(mp_state_seq(lf8_dtrans_cond_mup) == 1) = [];
%     lf8_dtrans_cond_mdown = lf8_dtrans;
%     lf8_dtrans_cond_mdown(mp_state_seq(lf8_dtrans_cond_mdown)==2) = [];
%     
    mp_utrig_mp(d,:) = get_event_trig_avg(wcv,mp_utrans,forwardlag,backwardlag);
    mp_dtrig_mp(d,:) = get_event_trig_avg(wcv,mp_dtrans,forwardlag,backwardlag);
    mp_utrig_lf8(d,:) = get_event_trig_avg(lf8,mp_utrans,forwardlag,backwardlag);
    mp_dtrig_lf8(d,:) = get_event_trig_avg(lf8,mp_dtrans,forwardlag,backwardlag);
    mp_utrig_spk(d,:) = get_event_trig_avg(bin_spk,mp_utrans,forwardlag,backwardlag);
    mp_dtrig_spk(d,:) = get_event_trig_avg(bin_spk,mp_dtrans,forwardlag,backwardlag);
    mp_utrig_lf2rlf(d,:) = get_event_trig_avg(lf2_rlf,mp_utrans,forwardlag,backwardlag);
    mp_dtrig_lf2rlf(d,:) = get_event_trig_avg(lf2_rlf,mp_dtrans,forwardlag,backwardlag);
    mp_utrig_lf2rhf(d,:) = get_event_trig_avg(lf2_rhf,mp_utrans,forwardlag,backwardlag);
    mp_dtrig_lf2rhf(d,:) = get_event_trig_avg(lf2_rhf,mp_dtrans,forwardlag,backwardlag);
    mp_utrig_lf2hf(d,:) = get_event_trig_avg(lf2_hf,mp_utrans,forwardlag,backwardlag);
    mp_dtrig_lf2hf(d,:) = get_event_trig_avg(lf2_hf,mp_dtrans,forwardlag,backwardlag);
    lf8_utrig_lf2rlf(d,:) = get_event_trig_avg(lf2_rlf,lf8_utrans,forwardlag,backwardlag);
    lf8_dtrig_lf2rlf(d,:) = get_event_trig_avg(lf2_rlf,lf8_dtrans,forwardlag,backwardlag);
    lf8_utrig_lf2rhf(d,:) = get_event_trig_avg(lf2_rhf,lf8_utrans,forwardlag,backwardlag);
    lf8_dtrig_lf2rhf(d,:) = get_event_trig_avg(lf2_rhf,lf8_dtrans,forwardlag,backwardlag);
    lf8_utrig_lf2hf(d,:) = get_event_trig_avg(lf2_hf,lf8_utrans,forwardlag,backwardlag);
    lf8_dtrig_lf2hf(d,:) = get_event_trig_avg(lf2_hf,lf8_dtrans,forwardlag,backwardlag);
    lf8_utrig_lf8(d,:) = get_event_trig_avg(lf8,lf8_utrans,forwardlag,backwardlag);
    lf8_dtrig_lf8(d,:) = get_event_trig_avg(lf8,lf8_dtrans,forwardlag,backwardlag);
    lf8_utrig_mp(d,:) = get_event_trig_avg(wcv,lf8_utrans,forwardlag,backwardlag);
    lf8_dtrig_mp(d,:) = get_event_trig_avg(wcv,lf8_dtrans,forwardlag,backwardlag);
    lf8_utrig_spk(d,:) = get_event_trig_avg(bin_spk,lf8_utrans,forwardlag,backwardlag);
    lf8_dtrig_spk(d,:) = get_event_trig_avg(bin_spk,lf8_dtrans,forwardlag,backwardlag);
%     lf8_utrig_lf2r_cond_mup(d,:) = get_event_trig_avg(lf2_r,lf8_utrans_cond_mup,forwardlag,backwardlag);
%     lf8_utrig_lf2r_cond_mdown(d,:) = get_event_trig_avg(lf2_r,lf8_utrans_cond_mdown,forwardlag,backwardlag);
%     lf8_dtrig_lf2r_cond_mup(d,:) = get_event_trig_avg(lf2_r,lf8_dtrans_cond_mup,forwardlag,backwardlag);
%     lf8_dtrig_lf2r_cond_mdown(d,:) = get_event_trig_avg(lf2_r,lf8_dtrans_cond_mdown,forwardlag,backwardlag);
%     lf8_utrig_mp_cond_mup(d,:) = get_event_trig_avg(wcv,lf8_utrans_cond_mup,forwardlag,backwardlag);
%     lf8_utrig_mp_cond_mdown(d,:) = get_event_trig_avg(wcv,lf8_utrans_cond_mdown,forwardlag,backwardlag);
%     lf8_dtrig_mp_cond_mup(d,:) = get_event_trig_avg(wcv,lf8_dtrans_cond_mup,forwardlag,backwardlag);
%     lf8_dtrig_mp_cond_mdown(d,:) = get_event_trig_avg(wcv,lf8_dtrans_cond_mdown,forwardlag,backwardlag);
%     lf8_utrig_spk_cond_mup(d,:) = get_event_trig_avg(bin_spk,lf8_utrans_cond_mup,forwardlag,backwardlag);
%     lf8_utrig_spk_cond_mdown(d,:) = get_event_trig_avg(bin_spk,lf8_utrans_cond_mdown,forwardlag,backwardlag);
%     lf8_dtrig_spk_cond_mup(d,:) = get_event_trig_avg(bin_spk,lf8_dtrans_cond_mup,forwardlag,backwardlag);
%     lf8_dtrig_spk_cond_mdown(d,:) = get_event_trig_avg(bin_spk,lf8_dtrans_cond_mdown,forwardlag,backwardlag);
    
end

%%
cd E:\WC_Germany\overall_EC\
save overall_EC_trig_lf2lfhf_layer3_data lags lf8_* mp_*

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
% 
% %%
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
