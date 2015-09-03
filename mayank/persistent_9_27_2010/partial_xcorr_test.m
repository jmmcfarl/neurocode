clear all
close all

load G:\WC_Germany\overall_EC\overall_EC_dir
addpath('G:\WC_Germany\parietal_cortical_2010\')
addpath('G:\WC_Germany\hsmm_state_detection\')

cd G:\WC_Germany\persistent_9_27_2010\
load pa_corresponding_lfp_revised_simp_new2

used_data = [l3mec_p l3lec_p];
sess_data = sess_data(used_data);

drive_letter = 'G';

dsf = 32;
Fsd = 2016/dsf;
minSegLength = 60;
maxLag = 2*Fsd;
niqf = 2016/2;
lcf = .05/niqf;
hcf = 10/niqf;
[b,a] = butter(2,[lcf hcf]);

mec = 1:22;
lec = 23:36;
bad_lf3 = [7 9 11 12 13 16 17 25 31 32 35];
good_lf3 = setdiff(1:36,bad_lf3);
good_lf3_mec = good_lf3(ismember(good_lf3,mec));
good_lf3_lec = good_lf3(ismember(good_lf3,lec));

% for i = 1:length(sess_data)
%     g4(i) = sess_data(i).gains(4);
% end
% sess_data(isnan(g4)) = [];

for d = 1:length(sess_data)
    %     d=14
    
    cdir = sess_data(d).directory;
    cdir(1) = 'G';
    disp(sprintf('session %d',d))
    cd(cdir);
    s_name = strcat(sess_data(d).region,'_l',sess_data(d).layer,'_',sess_data(d).name);
    
    load ./desynch_times_lf8
    load ./used_data wcv_minus_spike lf8 lf3 lf5 lf2
    lf3_g = lf3 + lf5;
    
%     lf3 = lf2;
    
    wcv_f = filtfilt(b,a,wcv_minus_spike);
    lf8_f = filtfilt(b,a,lf8);
%     lf3_f = filtfilt(b,a,lf3);
%     lf2_f = filtfilt(b,a,lf2);
%     lf4_f = filtfilt(b,a,lf4);
    
    down_w = downsample(wcv_f,dsf);
    down_8 = downsample(lf8_f,dsf);
%     down_3 = downsample(lf3_f,dsf);
%     down_2 = downsample(lf2_f,dsf);
%     down_4 = downsample(lf4_f,dsf);
    
%     csd = 2*down_3 - down_2 - down_4;
    
    down_3 = get_hf_features(lf3_g,2016,Fsd,[15 80],0.05);
%     down_3 = get_hf_features(lf3,2016,Fsd,[10 80],0.075);
    
    %zscore
    down_w = zscore(down_w);
    down_8 = zscore(down_8);
    down_3 = zscore(down_3);
%     down_2 = zscore(down_2);
%     csd = zscore(csd);
    
    
    load ./pa_hsmm_state_seq8_new2
    load ./pa_hsmm_state_seq_new2
    lfp_state_seq = hsmm_bbstate_seq8;
    mp_state_seq = hsmm_bbstate_seq;
    
    [new_seg_inds] = resample_uds_seg_inds(hmm8.UDS_segs,hmm8.Fs,252,2*length(down_8));
    [up_trans_inds8,down_trans_inds8] = compute_state_transitions_seg(new_seg_inds,lfp_state_seq);
    [up_trans_inds,down_trans_inds] = compute_state_transitions_seg(new_seg_inds,mp_state_seq);
    up_trans_inds8 = round(up_trans_inds8/4);
    up_trans_inds = round(up_trans_inds/4);
    down_trans_inds8 = round(down_trans_inds8/4);
    down_trans_inds = round(down_trans_inds/4);
    
    %     skipped_lfp_trans = [];
    %     for j = 1:length(up_trans_inds)
    %         skipped_lfp_trans = [skipped_lfp_trans mp_upskipped{d}.inds{j}];
    %     end
    %     non_skipped_lfp_trans = setdiff(1:length(down_trans_inds8),skipped_lfp_trans);
    %     n_uskip(d) = length(skipped_lfp_trans);
    %
    %     dskipped_lfp_trans = [];
    %     for j = 1:length(up_trans_inds)
    %         dskipped_lfp_trans = [dskipped_lfp_trans mp_downskipped{d}.inds{j}];
    %     end
    %     non_dskipped_lfp_trans = setdiff(1:length(up_trans_inds8),dskipped_lfp_trans);
    %     n_dskip(d) = length(dskipped_lfp_trans);
    
    up_trans_vec8 = zeros(size(down_8));
    up_trans_vec8(up_trans_inds8) = 1;
    %     up_trans_vec = zeros(size(down_8));
    %     up_trans_vec(up_trans_inds) = 1;
    %     up_trans_vec8_ns = zeros(size(down_8));
    %     up_trans_vec8_ns(up_trans_inds8(non_skipped_lfp_trans)) = 1;
    %     up_trans_vec8_s = zeros(size(down_8));
    %     up_trans_vec8_s(up_trans_inds8(skipped_lfp_trans)) = 1;
    up_trans_vec = zeros(size(down_8));
    up_trans_vec(up_trans_inds) = 1;
    %     down_trans_vec8_ns = zeros(size(down_8));
    %     down_trans_vec8_ns(down_trans_inds8(non_dskipped_lfp_trans)) = 1;
    %     down_trans_vec8_s = zeros(size(down_8));
    %     down_trans_vec8_s(down_trans_inds8(dskipped_lfp_trans)) = 1;
    down_trans_vec8 = zeros(size(down_8));
    down_trans_vec8(down_trans_inds8) = 1;
    down_trans_vec = zeros(size(down_8));
    down_trans_vec(down_trans_inds) = 1;
    
    Xwd = zeros(2*maxLag+1,length(down_8));
    Xwu = zeros(2*maxLag+1,length(down_8));
    for i = 1:maxLag
        Xwd(i,1:end-maxLag+i) = down_trans_vec(maxLag-i+1:end);
        Xwu(i,1:end-maxLag+i) = up_trans_vec(maxLag-i+1:end);
    end
    for i = (maxLag+1):(2*maxLag+1)
        Xwd(i,(i-maxLag+1):end) = down_trans_vec(1:end-i+maxLag);
        Xwu(i,(i-maxLag+1):end) = up_trans_vec(1:end-i+maxLag);
    end
    X8d = zeros(2*maxLag+1,length(down_8));
    X8u = zeros(2*maxLag+1,length(down_8));
    for i = 1:maxLag
        X8d(i,1:end-maxLag+i) = down_trans_vec8(maxLag-i+1:end);
        X8u(i,1:end-maxLag+i) = up_trans_vec8(maxLag-i+1:end);
    end
    for i = (maxLag+1):(2*maxLag+1)
        X8d(i,(i-maxLag+1):end) = down_trans_vec8(1:end-i+maxLag);
        X8u(i,(i-maxLag+1):end) = up_trans_vec8(1:end-i+maxLag);
    end
    
    %     X8d_ns = zeros(2*maxLag+1,length(down_8));
    %     X8d_s = zeros(2*maxLag+1,length(down_8));
    %     X8u = zeros(2*maxLag+1,length(down_8));
    %     for i = 1:maxLag
    %         X8d_ns(i,1:end-maxLag+i) = down_trans_vec8_ns(maxLag-i+1:end);
    %         X8d_s(i,1:end-maxLag+i) = down_trans_vec8_s(maxLag-i+1:end);
    %         X8u(i,1:end-maxLag+i) = up_trans_vec8_ns(maxLag-i+1:end);
    %     end
    %     for i = (maxLag+1):(2*maxLag+1)
    %         X8d_ns(i,(i-maxLag+1):end) = down_trans_vec8_ns(1:end-i+maxLag);
    %         X8d_s(i,(i-maxLag+1):end) = down_trans_vec8_s(1:end-i+maxLag);
    %         X8u(i,(i-maxLag+1):end) = up_trans_vec8_ns(1:end-i+maxLag);
    %     end
    
    %
    %compute markers indicating segments of data to be used
    t_axis = (1:length(down_8))/Fsd;
    if ~isempty(desynch_times_lf8)
        desynch_start = round(interp1(t_axis,1:length(t_axis),desynch_times_lf8(:,1)));
        desynch_stop = round(interp1(t_axis,1:length(t_axis),desynch_times_lf8(:,2)));
    else
        desynch_start = [];
        desynch_stop = [];
    end
    desynch_ind = zeros(size(down_8));
    for i = 1:length(desynch_start)
        desynch_ind(desynch_start(i):desynch_stop(i)) = 1;
    end
    synch_starts = find(desynch_ind(1:end-1)==1 & desynch_ind(2:end)==0)+1;
    if desynch_ind(1) == 0
        synch_starts = [1; synch_starts];
    end
    synch_stops = find(desynch_ind(1:end-1)==0 & desynch_ind(2:end)==1)+1;
    if desynch_ind(end) == 0
        synch_stops = [synch_stops; length(down_8)];
    end
    sMarkers = [synch_starts(:) synch_stops(:)];
    
    seg_durs = diff(sMarkers')'/Fsd;
    too_short = find(seg_durs < minSegLength);
    sMarkers(too_short,:) = [];
    seg_durs(too_short) = [];
    %
    cnt = 0;
    for i = 1:size(sMarkers,1)
        cnt = cnt+1;
        
%         B(cnt,:) = regress(down_3(sMarkers(i,1):sMarkers(i,2)),...
%             [Xw(:,sMarkers(i,1):sMarkers(i,2))' X8(:,sMarkers(i,1):sMarkers(i,2))']);
%         Bw(cnt,:) = regress(down_3(sMarkers(i,1):sMarkers(i,2)),...
%             [Xw(:,sMarkers(i,1):sMarkers(i,2))']);
%         B8(cnt,:) = regress(down_3(sMarkers(i,1):sMarkers(i,2)),...
%             [X8(:,sMarkers(i,1):sMarkers(i,2))']);
        % cL = length(sMarkers(i,1):sMarkers(i,2));
        % [B(cnt,:),~,~,~,stats] = regress(down_3(sMarkers(i,1):sMarkers(i,2)),...
        %     [ones(1,cL)' Xwd(:,sMarkers(i,1):sMarkers(i,2))' X8d(:,sMarkers(i,1):sMarkers(i,2))' ...
        %     Xwu(:,sMarkers(i,1):sMarkers(i,2))' X8u(:,sMarkers(i,1):sMarkers(i,2))']);
        % R2(cnt) = stats(1);
        % rerr(cnt) = stats(4);
        % [B_mwd(cnt,:),~,~,~,stats] = regress(down_3(sMarkers(i,1):sMarkers(i,2)),...
        %     [ones(1,cL)' X8d(:,sMarkers(i,1):sMarkers(i,2))' ...
        %     Xwu(:,sMarkers(i,1):sMarkers(i,2))' X8u(:,sMarkers(i,1):sMarkers(i,2))']);
        % R2_mwd(cnt) = stats(1);
        % rerr_mwd(cnt) = stats(4);
        % [B_mwu(cnt,:),~,~,~,stats] = regress(down_3(sMarkers(i,1):sMarkers(i,2)),...
        %     [ones(1,cL)' Xwd(:,sMarkers(i,1):sMarkers(i,2))' X8d(:,sMarkers(i,1):sMarkers(i,2))' ...
        %     X8u(:,sMarkers(i,1):sMarkers(i,2))']);
        % R2_mwu(cnt) = stats(1);
        % rerr_mwu(cnt) = stats(4);
        
        
        [Bw(cnt,:),~,~,~,statsw(cnt,:)] = regress(down_3(sMarkers(i,1):sMarkers(i,2)),...
            [Xwd(:,sMarkers(i,1):sMarkers(i,2))' Xwu(:,sMarkers(i,1):sMarkers(i,2))']);
        [B8(cnt,:),~,~,~,stats8(cnt,:)] = regress(down_3(sMarkers(i,1):sMarkers(i,2)),...
            [X8d(:,sMarkers(i,1):sMarkers(i,2))' X8u(:,sMarkers(i,1):sMarkers(i,2))']);
        
%         B(cnt,:) = regress(down_3(sMarkers(i,1):sMarkers(i,2)),...
%             [Xwd(:,sMarkers(i,1):sMarkers(i,2))' Xwu(:,sMarkers(i,1):sMarkers(i,2))' ...
%             X8u(:,sMarkers(i,1):sMarkers(i,2))' ...
%             X8d_ns(:,sMarkers(i,1):sMarkers(i,2))' X8d_s(:,sMarkers(i,1):sMarkers(i,2))']);
        B(cnt,:) = regress(down_3(sMarkers(i,1):sMarkers(i,2)),...
            [Xwd(:,sMarkers(i,1):sMarkers(i,2))' Xwu(:,sMarkers(i,1):sMarkers(i,2))' ...
            X8u(:,sMarkers(i,1):sMarkers(i,2))' X8d(:,sMarkers(i,1):sMarkers(i,2))']);
%         B8(cnt,:) = regress(down_3(sMarkers(i,1):sMarkers(i,2)),...
%             [X8u(:,sMarkers(i,1):sMarkers(i,2))' X8d_ns(:,sMarkers(i,1):sMarkers(i,2))' ...
%             X8d_s(:,sMarkers(i,1):sMarkers(i,2))']);
        
    end
    %
    lags = -maxLag:maxLag;
    NLags = length(lags);
    
    %compute weighted averages (weighted by relative duration
    weighting = seg_durs/sum(seg_durs);

    B_avg = sum(B.*repmat(weighting(:),1,size(B,2)),1);
%     Bmwd_avg = sum(B_mwd.*repmat(weighting(:),1,size(B_mwd,2)),1);
%     Bmwu_avg = sum(B_mwu.*repmat(weighting(:),1,size(B_mwu,2)),1);
    %     R2_avg(d) = sum(R2(:).*weighting(:));
    %     R2_mwd_avg(d) = sum(R2_mwd(:).*weighting(:));
    %     R2_mwu_avg(d) = sum(R2_mwu(:).*weighting(:));
    %     rerr_avg(d) = sum(rerr(:).*weighting(:));
    %     rerr_mwd_avg(d) = sum(rerr_mwd(:).*weighting(:));
    %     rerr_mwu_avg(d) = sum(rerr_mwu(:).*weighting(:));
    
    Bw_avg = sum(Bw.*repmat(weighting(:),1,size(Bw,2)),1);
    B8_avg = sum(B8.*repmat(weighting(:),1,size(B8,2)),1);
    
%     C3w_part(d,:) = B_avg(1:length(lags));
%     C38_part(d,:) = B_avg(length(lags)+1:end);
%     C3w(d,:) = Bw_avg;
%     C38(d,:) = B8_avg;
    
    C3w_d(d,:) = Bw_avg(1:NLags);
    C3w_u(d,:) = Bw_avg(NLags+1:end);
    C38_d(d,:) = B8_avg(1:NLags);
    C38_u(d,:) = B8_avg(NLags+1:end);
%     C3w_dpart(d,:) = B_avg(1:NLags);
%     C38_dpart(d,:) = B_avg(NLags+1:2*NLags);
%     C3w_upart(d,:) = B_avg(2*NLags+1:3*NLags);
%     C38_upart(d,:) = B_avg(3*NLags+1:end);
    
%     C38_u(d,:) = B8_avg(1:NLags);
%     C38_d_ns(d,:) = B8_avg(NLags+1:2*NLags);
%     C38_d_s(d,:) = B8_avg(2*NLags+1:end);
    C3w_dpart(d,:) = B_avg(1:NLags);
    C3w_upart(d,:) = B_avg(NLags+1:2*NLags);
    C38_upart(d,:) = B_avg(2*NLags+1:3*NLags);
    C38_dpart(d,:) = B_avg(3*NLags+1:4*NLags);
%     C38_d_nspart(d,:) = B_avg(3*NLags+1:4*NLags);
%     C38_d_spart(d,:) = B_avg(4*NLags+1:end);
    
    clear down_w down_8 wcv* lf8* lf3* B* Xw X8 stats* R2 R2_mwd R2_mwu rerr rerr_mwd rerr_mwu
    
end
%%

cd G:\WC_Germany\persistent_9_27_2010
save partial_sta_analysis_d24_fulltrans_lf3ghf lags Fsd C*

%%
figure
h = errorbar(lags/Fsd,nanmean(C38_d_ns(good_lf3_mec,:)),nanstd(C38_d_ns(good_lf3_mec,:))/sqrt(length(good_lf3_mec)));
errorbar_tick(h,.001,'units');
hold on
h = errorbar(lags/Fsd,nanmean(C38_d_s(good_lf3_mec,:)),nanstd(C38_d_s(good_lf3_mec,:))/sqrt(length(good_lf3_mec)),'r');
errorbar_tick(h,.001,'units');
h = errorbar(lags/Fsd,nanmean(C38_d_nspart(good_lf3_mec,:)),nanstd(C38_d_nspart(good_lf3_mec,:))/sqrt(length(good_lf3_mec)),'k');
errorbar_tick(h,.001,'units');
h = errorbar(lags/Fsd,nanmean(C38_d_spart(good_lf3_mec,:)),nanstd(C38_d_spart(good_lf3_mec,:))/sqrt(length(good_lf3_mec)),'g');
errorbar_tick(h,.001,'units');
xlim([-2 2])

%%
figure
h = errorbar(lags/Fsd,nanmean(C3w(good_lf3_mec,:)),nanstd(C3w(good_lf3_mec,:))/sqrt(length(good_lf3_mec)));
errorbar_tick(h,.001,'units');
hold on
h = errorbar(lags/Fsd,nanmean(C38(good_lf3_mec,:)),nanstd(C38(good_lf3_mec,:))/sqrt(length(good_lf3_mec)),'r');
errorbar_tick(h,.001,'units');
h = errorbar(lags/Fsd,nanmean(C3w_part(good_lf3_mec,:)),nanstd(C3w_part(good_lf3_mec,:))/sqrt(length(good_lf3_mec)),'g');
errorbar_tick(h,.001,'units');
h = errorbar(lags/Fsd,nanmean(C38_part(good_lf3_mec,:)),nanstd(C38_part(good_lf3_mec,:))/sqrt(length(good_lf3_mec)),'g');
errorbar_tick(h,.001,'units');

figure
h = errorbar(lags/Fsd,nanmean(C3w(good_lf3_lec,:)),nanstd(C3w(good_lf3_lec,:))/sqrt(length(good_lf3_lec)));
errorbar_tick(h,.001,'units');
hold on
h = errorbar(lags/Fsd,nanmean(C38(good_lf3_lec,:)),nanstd(C38(good_lf3_lec,:))/sqrt(length(good_lf3_lec)),'r');
errorbar_tick(h,.001,'units');
h = errorbar(lags/Fsd,nanmean(C3w_part(good_lf3_lec,:)),nanstd(C3w_part(good_lf3_lec,:))/sqrt(length(good_lf3_lec)),'g');
errorbar_tick(h,.001,'units');
h = errorbar(lags/Fsd,nanmean(C38_part(good_lf3_lec,:)),nanstd(C38_part(good_lf3_lec,:))/sqrt(length(good_lf3_lec)),'g');
errorbar_tick(h,.001,'units');


%%
figure
h = errorbar(lags/Fsd,nanmean(C3w_u(good_lf3_mec,:)),nanstd(C3w_u(good_lf3_mec,:))/sqrt(length(good_lf3_mec)),'b');
errorbar_tick(h,.001,'units');
hold on
h = errorbar(lags/Fsd,nanmean(C38_u(good_lf3_mec,:)),nanstd(C38_u(good_lf3_mec,:))/sqrt(length(good_lf3_mec)),'r');
errorbar_tick(h,.001,'units');
h = errorbar(lags/Fsd,nanmean(C3w_upart(good_lf3_mec,:)),nanstd(C3w_upart(good_lf3_mec,:))/sqrt(length(good_lf3_mec)),'k');
errorbar_tick(h,.001,'units');
h = errorbar(lags/Fsd,nanmean(C38_upart(good_lf3_mec,:)),nanstd(C38_upart(good_lf3_mec,:))/sqrt(length(good_lf3_mec)),'c');
errorbar_tick(h,.001,'units');
xlim([-2 2])

figure
h = errorbar(lags/Fsd,nanmean(C3w_u(good_lf3_lec,:)),nanstd(C3w_u(good_lf3_lec,:))/sqrt(length(good_lf3_lec)),'b');
errorbar_tick(h,.001,'units');
hold on
h = errorbar(lags/Fsd,nanmean(C38_u(good_lf3_lec,:)),nanstd(C38_u(good_lf3_lec,:))/sqrt(length(good_lf3_lec)),'r');
errorbar_tick(h,.001,'units');
h = errorbar(lags/Fsd,nanmean(C3w_upart(good_lf3_lec,:)),nanstd(C3w_upart(good_lf3_lec,:))/sqrt(length(good_lf3_lec)),'k');
errorbar_tick(h,.001,'units');
h = errorbar(lags/Fsd,nanmean(C38_upart(good_lf3_lec,:)),nanstd(C38_upart(good_lf3_lec,:))/sqrt(length(good_lf3_lec)),'c');
errorbar_tick(h,.001,'units');
xlim([-2 2])

figure
h = errorbar(lags/Fsd,nanmean(C3w_d(good_lf3_mec,:)),nanstd(C3w_d(good_lf3_mec,:))/sqrt(length(good_lf3_mec)),'b');
errorbar_tick(h,.001,'units');
hold on
h = errorbar(lags/Fsd,nanmean(C38_d(good_lf3_mec,:)),nanstd(C38_u(good_lf3_mec,:))/sqrt(length(good_lf3_mec)),'r');
errorbar_tick(h,.001,'units');
h = errorbar(lags/Fsd,nanmean(C3w_dpart(good_lf3_mec,:)),nanstd(C3w_dpart(good_lf3_mec,:))/sqrt(length(good_lf3_mec)),'k');
errorbar_tick(h,.001,'units');
h = errorbar(lags/Fsd,nanmean(C38_dpart(good_lf3_mec,:)),nanstd(C38_dpart(good_lf3_mec,:))/sqrt(length(good_lf3_mec)),'c');
errorbar_tick(h,.001,'units');
xlim([-2 2])

figure
h = errorbar(lags/Fsd,nanmean(C3w_d(good_lf3_lec,:)),nanstd(C3w_d(good_lf3_lec,:))/sqrt(length(good_lf3_lec)),'c');
errorbar_tick(h,.001,'units');
hold on
h = errorbar(lags/Fsd,nanmean(C38_d(good_lf3_lec,:)),nanstd(C38_u(good_lf3_lec,:))/sqrt(length(good_lf3_lec)),'c');
errorbar_tick(h,.001,'units');
h = errorbar(lags/Fsd,nanmean(C3w_dpart(good_lf3_lec,:)),nanstd(C3w_dpart(good_lf3_lec,:))/sqrt(length(good_lf3_lec)),'k');
errorbar_tick(h,.001,'units');
h = errorbar(lags/Fsd,nanmean(C38_dpart(good_lf3_lec,:)),nanstd(C38_dpart(good_lf3_lec,:))/sqrt(length(good_lf3_lec)),'k');
errorbar_tick(h,.001,'units');
xlim([-2 2])

