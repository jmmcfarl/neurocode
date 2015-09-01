clear all
close all

load G:\WC_Germany\overall_EC\overall_EC_dir
addpath('G:\WC_Germany\parietal_cortical_2010\')
addpath('G:\WC_Germany\hsmm_state_detection\')

used_data = [l3mec_p l3lec_p];
sess_data = sess_data(used_data);

drive_letter = 'G';

dsf = 8;
Fsd = 2016/dsf;
minSegLength = 60;
maxLag = 10*Fsd;
niqf = 2016/2;
lcf = .05/niqf;
hcf = 2/niqf;
[b,a] = butter(2,[lcf hcf]);

for d = 1:length(sess_data)

    cdir = sess_data(d).directory;
    cdir(1) = 'G';
    disp(sprintf('session %d',d))
    cd(cdir);
    s_name = strcat(sess_data(d).region,'_l',sess_data(d).layer,'_',sess_data(d).name);
    
    load ./desynch_times_lf8
    load ./used_data wcv_minus_spike lf8 lf3 lf5 lf2

    lf2_r = lf2 - lf5;
    wcv_f = filtfilt(b,a,wcv_minus_spike);
    lf8_f = filtfilt(b,a,lf8);
    lf5_f = filtfilt(b,a,lf5);
    lf3_f = filtfilt(b,a,lf3);
    lf2_f = filtfilt(b,a,lf2);
    lf2_rf = filtfilt(b,a,lf2_r);
    
    down_w = downsample(wcv_f,dsf);
    down_8 = downsample(lf8_f,dsf);
    down_5 = downsample(lf5_f,dsf);
    down_3 = downsample(lf3_f,dsf);
    down_2 = downsample(lf2_f,dsf);
    down_2r = downsample(lf2_rf,dsf);
    
    %zscore
    down_w = zscore(down_w);
    down_8 = zscore(down_8);
    down_5 = zscore(down_5);
    down_3 = zscore(down_3);
    down_2 = zscore(down_2);
    down_2r = zscore(down_2r);
    
    down_2h = get_hf_features(lf2,2016,Fsd,[20 80],.05);
    down_3h = get_hf_features(lf3,2016,Fsd,[20 80],.05);
    down_5h = get_hf_features(lf5,2016,Fsd,[20 80],.05);
    down_8h = get_hf_features(lf8,2016,Fsd,[20 80],.05);
    
    
    %compute markers indicating segments of data to be used
    t_axis = (1:length(down_w))/Fsd;
    if ~isempty(desynch_times_lf8)
        desynch_start = round(interp1(t_axis,1:length(t_axis),desynch_times_lf8(:,1)));
        desynch_stop = round(interp1(t_axis,1:length(t_axis),desynch_times_lf8(:,2)));
    else
        desynch_start = [];
        desynch_stop = [];
    end
    desynch_ind = zeros(size(down_w));
    for i = 1:length(desynch_start)
        desynch_ind(desynch_start(i):desynch_stop(i)) = 1;
    end
    synch_starts = find(desynch_ind(1:end-1)==1 & desynch_ind(2:end)==0)+1;
    if desynch_ind(1) == 0
        synch_starts = [1; synch_starts];
    end
    synch_stops = find(desynch_ind(1:end-1)==0 & desynch_ind(2:end)==1)+1;
    if desynch_ind(end) == 0
        synch_stops = [synch_stops; length(down_w)];
    end
    sMarkers = [synch_starts(:) synch_stops(:)];
    
    seg_durs = diff(sMarkers')'/Fsd;
    too_short = find(seg_durs < minSegLength);
    sMarkers(too_short,:) = [];
    seg_durs(too_short) = [];
    
    cnt = 0;
    for i = 1:size(sMarkers,1)
        cnt = cnt+1;
        
        [wcv_acorr(cnt,:),lags] = xcov(down_w(sMarkers(i,1):sMarkers(i,2)),maxLag,'coeff');
        [lf8_acorr(cnt,:),lags] = xcov(down_8(sMarkers(i,1):sMarkers(i,2)),maxLag,'coeff');
        [lf3_acorr(cnt,:),lags] = xcov(down_3(sMarkers(i,1):sMarkers(i,2)),maxLag,'coeff');
        
        [w8_x(cnt,:),lags] = xcov(down_w(sMarkers(i,1):sMarkers(i,2)),down_8(sMarkers(i,1):sMarkers(i,2)),maxLag,'coeff');
        [x53(cnt,:),lags] = xcov(down_5(sMarkers(i,1):sMarkers(i,2)),down_3(sMarkers(i,1):sMarkers(i,2)),maxLag,'coeff');
        [x83(cnt,:),lags] = xcov(down_8(sMarkers(i,1):sMarkers(i,2)),down_3(sMarkers(i,1):sMarkers(i,2)),maxLag,'coeff');
        [x82(cnt,:),lags] = xcov(down_8(sMarkers(i,1):sMarkers(i,2)),down_2(sMarkers(i,1):sMarkers(i,2)),maxLag,'coeff');
        [x82r(cnt,:),lags] = xcov(down_8(sMarkers(i,1):sMarkers(i,2)),down_2r(sMarkers(i,1):sMarkers(i,2)),maxLag,'coeff');
        [w3_x(cnt,:),lags] = xcov(down_w(sMarkers(i,1):sMarkers(i,2)),down_3(sMarkers(i,1):sMarkers(i,2)),maxLag,'coeff');
        [w2_x(cnt,:),lags] = xcov(down_w(sMarkers(i,1):sMarkers(i,2)),down_2(sMarkers(i,1):sMarkers(i,2)),maxLag,'coeff');

        [x88h(cnt,:),lags] = xcov(down_8(sMarkers(i,1):sMarkers(i,2)),down_8h(sMarkers(i,1):sMarkers(i,2)),maxLag,'coeff');
        [x85h(cnt,:),lags] = xcov(down_8(sMarkers(i,1):sMarkers(i,2)),down_5h(sMarkers(i,1):sMarkers(i,2)),maxLag,'coeff');
        [x83h(cnt,:),lags] = xcov(down_8(sMarkers(i,1):sMarkers(i,2)),down_3h(sMarkers(i,1):sMarkers(i,2)),maxLag,'coeff');
        [x82h(cnt,:),lags] = xcov(down_8(sMarkers(i,1):sMarkers(i,2)),down_2h(sMarkers(i,1):sMarkers(i,2)),maxLag,'coeff');

        [xw8h(cnt,:),lags] = xcov(down_w(sMarkers(i,1):sMarkers(i,2)),down_8h(sMarkers(i,1):sMarkers(i,2)),maxLag,'coeff');
        [xw5h(cnt,:),lags] = xcov(down_w(sMarkers(i,1):sMarkers(i,2)),down_5h(sMarkers(i,1):sMarkers(i,2)),maxLag,'coeff');
        [xw3h(cnt,:),lags] = xcov(down_w(sMarkers(i,1):sMarkers(i,2)),down_3h(sMarkers(i,1):sMarkers(i,2)),maxLag,'coeff');
        [xw2h(cnt,:),lags] = xcov(down_w(sMarkers(i,1):sMarkers(i,2)),down_2h(sMarkers(i,1):sMarkers(i,2)),maxLag,'coeff');

    end

    %compute weighted averages (weighted by relative duration
    weighting = seg_durs/sum(seg_durs);
    tot_wcv_acorr(d,:) = sum(wcv_acorr.*repmat(weighting(:),1,length(lags)),1);
    tot_lf8_acorr(d,:) = sum(lf8_acorr.*repmat(weighting(:),1,length(lags)),1);
    tot_lf3_acorr(d,:) = sum(lf3_acorr.*repmat(weighting(:),1,length(lags)),1);
    tot_w8_x(d,:) = sum(w8_x.*repmat(weighting(:),1,length(lags)),1);
    tot_53_x(d,:) = sum(x53.*repmat(weighting(:),1,length(lags)),1);
    tot_83_x(d,:) = sum(x83.*repmat(weighting(:),1,length(lags)),1);
    tot_82_x(d,:) = sum(x82.*repmat(weighting(:),1,length(lags)),1);
    tot_82r_x(d,:) = sum(x82r.*repmat(weighting(:),1,length(lags)),1);
    tot_w3_x(d,:) = sum(w3_x.*repmat(weighting(:),1,length(lags)),1);
    tot_w2_x(d,:) = sum(w2_x.*repmat(weighting(:),1,length(lags)),1);

    tot_88h_x(d,:) = sum(x88h.*repmat(weighting(:),1,length(lags)),1);
    tot_85h_x(d,:) = sum(x85h.*reptime_domain_data_lf2mat(weighting(:),1,length(lags)),1);
    tot_83h_x(d,:) = sum(x83h.*repmat(weighting(:),1,length(lags)),1);
    tot_82h_x(d,:) = sum(x82h.*repmat(weighting(:),1,length(lags)),1);
    
    tot_w8h_x(d,:) = sum(xw8h.*repmat(weighting(:),1,length(lags)),1);
    tot_w5h_x(d,:) = sum(xw5h.*repmat(weighting(:),1,length(lags)),1);
    tot_w3h_x(d,:) = sum(xw3h.*repmat(weighting(:),1,length(lags)),1);
    tot_w2h_x(d,:) = sum(xw2h.*repmat(weighting(:),1,length(lags)),1);
    
        clear down_w down_8 wcv* lf8* lf3* w8_x* w3_x* w2_x* w3s_x* x53 x83 x82* x*h

end

cd G:\WC_Germany\persistent_9_27_2010\
save time_domain_data_lf2 tot* Fsd lags

%%
mec = 1:22;
lec = 23:36;
bad_lf3 = [7 9 11 12 13 16 17 25 31 32 35];
good_lf3 = setdiff(1:36,bad_lf3);
good_lf3_mec = good_lf3(ismember(good_lf3,mec));
good_lf3_lec = good_lf3(ismember(good_lf3,lec));

%%
for i = 1:36
   plot(lags/Fsd,tot_w8_x(i,:),'b')
   hold on
   plot(lags/Fsd,tot_82r_x(i,:),'k')
   plot(lags/Fsd,tot_83_x(i,:),'r')
   plot(lags/Fsd,tot_82h_x(i,:),'g')
   plot(lags/Fsd,tot_83h_x(i,:),'c')
   plot(lags/Fsd,tot_85h_x(i,:),'y')
   xlim([-1.5 1.5])
   i
   pause
   clf
   
end

%%
h = errorbar(lags/Fsd,fliplr(mean(tot_w8_x(mec,:))),fliplr(std(tot_w8_x(mec,:)))/sqrt(length(mec)));
errorbar_tick(h,.001,'units');
xlim([-2 2])
hold on
h = errorbar(lags/Fsd,mean(tot_83_x(good_lf3,:)),std(tot_83_x(good_lf3,:))/sqrt(length(good_lf3)),'k');
errorbar_tick(h,.001,'units');
h = errorbar(lags/Fsd,mean(tot_82_x(good_lf3,:)),std(tot_82_x(good_lf3,:))/sqrt(length(good_lf3)),'color',[0.2 0.8 0.1]);
errorbar_tick(h,.001,'units');
h = errorbar(lags/Fsd,mean(tot_83h_x(good_lf3,:)),std(tot_83h_x(good_lf3,:))/sqrt(length(good_lf3)),'g');
errorbar_tick(h,.001,'units');
h = errorbar(lags/Fsd,mean(tot_85h_x(good_lf3,:)),std(tot_85h_x(good_lf3,:))/sqrt(length(good_lf3)),'c');
errorbar_tick(h,.001,'units');
h = errorbar(lags/Fsd,mean(tot_88h_x(good_lf3,:)),std(tot_88h_x(good_lf3,:))/sqrt(length(good_lf3)),'color',[0.8 0.2 0.1]);
errorbar_tick(h,.001,'units');

figure
h = errorbar(lags/Fsd,mean(tot_w8_x(mec,:)),std(tot_w8_x(mec,:))/sqrt(length(mec)),'r');
errorbar_tick(h,.001,'units');
xlim([-2 2])
hold on
h = errorbar(lags/Fsd,mean(tot_w3_x(good_lf3_mec,:)),std(tot_w3_x(good_lf3_mec,:))/sqrt(length(good_lf3_mec)),'k');
errorbar_tick(h,.001,'units');
h = errorbar(lags/Fsd,mean(tot_w2_x(good_lf3_mec,:)),std(tot_w2_x(good_lf3_mec,:))/sqrt(length(good_lf3_mec)),'color',[0.2 0.8 0.1]);
errorbar_tick(h,.001,'units');
h = errorbar(lags/Fsd,mean(tot_w3h_x(good_lf3_mec,:)),std(tot_83h_x(good_lf3_mec,:))/sqrt(length(good_lf3_mec)),'g');
errorbar_tick(h,.001,'units');
h = errorbar(lags/Fsd,mean(tot_w5h_x(good_lf3_mec,:)),std(tot_w5h_x(good_lf3_mec,:))/sqrt(length(good_lf3_mec)),'c');
errorbar_tick(h,.001,'units');
h = errorbar(lags/Fsd,mean(tot_w8h_x(good_lf3_mec,:)),std(tot_w8h_x(good_lf3_mec,:))/sqrt(length(good_lf3_mec)),'color',[0.8 0.2 0.1]);
errorbar_tick(h,.001,'units');

yl = ylim();
line([0 0],yl,'color','k')