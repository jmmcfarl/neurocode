clear all
close all
cd C:\WC_Germany\Cortical_analysis

load C:\WC_Germany\Cortical_analysis\cortical_dir
load C:\WC_Germany\Cortical_analysis\desynch_detect\desynch_times
load C:\WC_Germany\Cortical_analysis\UDS_dur_raw\UDS_raw_data
load C:\WC_Germany\Cortical_analysis\UDS_synch_state_dur\UDS_synch_state_dur_data

dsf = 8;
Fsd = 2016/dsf;
minSegLength = 100;
maxLag = 10*Fsd;
niqf = 2016/2;
lcf = .05/niqf;
hcf = 40/niqf;
hcf2 = 5/niqf;

[b,a] = butter(2,[lcf hcf]);
[b2,a2] = butter(2,[lcf hcf2]);
for d = 1:length(sess_data)

    disp(sprintf('session %d',d))
    cd(sess_data(d).directory);
    pwd

    load used_data wcv_minus_spike lf8 

    wcv_f = filtfilt(b,a,wcv_minus_spike);
    lf8_f = filtfilt(b,a,lf8);
%     lf5_f = filtfilt(b,a,lf5);
     wcv_f2 = filtfilt(b2,a2,wcv_minus_spike);
    lf8_f2 = filtfilt(b2,a2,lf8);
   
    down_w = downsample(wcv_f,dsf);
    down_8 = downsample(lf8_f,dsf);
        down_w2 = downsample(wcv_f2,dsf);
    down_82 = downsample(lf8_f2,dsf);

%      down_5 = downsample(lf5_f,dsf);
   
    %zscor
    deriv_w = [0;diff(down_w2)];
    deriv_8 = [0;diff(down_82)];
    
    
    down_w = zscore(down_w);
    down_8 = zscore(down_8);
%     down_5 = zscore(down_5);
    
    deriv_w = zscore(deriv_w);
    deriv_8 = zscore(deriv_8);
    deriv_w(deriv_w < 0) = 0;
    deriv_8(deriv_8 < 0) = 0;
    
    desynch_start = desynch_start_times{d}*Fsd;
    desynch_stop = desynch_stop_times{d}*Fsd;

    if ~isempty(desynch_start)
        sMarkers = [[1;desynch_stop'] [desynch_start';length(down_w)]]
    else
        sMarkers = [1 length(down_w)];
    end

    cnt = 0;
    for i = 1:size(sMarkers,1)
        if (sMarkers(i,2)-sMarkers(i,1))/Fsd > minSegLength
            cnt = cnt+1;
            
            [wcv_acorr(cnt,:),lags] = xcov(down_w,maxLag,'coeff');
            [lf8_acorr(cnt,:),lags] = xcov(down_8,maxLag,'coeff');
            [wcv_acorr_d(cnt,:),lags] = xcov(deriv_w,maxLag,'coeff');
            [lf8_acorr_d(cnt,:),lags] = xcov(deriv_8,maxLag,'coeff');

%             [lf5_acorr(cnt,:),lags] = xcov(down_5,maxLag,'coeff');
            
            [w8_x(cnt,:),lags] = xcov(down_w,down_8,maxLag,'coeff');
            [w8_x_d(cnt,:),lags] = xcov(deriv_w,deriv_8,maxLag,'coeff');

%             [w5_x(cnt,:),lags] = xcov(down_w,down_5,maxLag,'coeff');
%             [l85_x(cnt,:),lags] = xcov(down_8,down_5,maxLag,'coeff');

        end
    end

    tot_wcv_acorr(d,:) = mean(wcv_acorr,1);
    tot_lf8_acorr(d,:) = mean(lf8_acorr,1);
    tot_wcv_acorr_d(d,:) = mean(wcv_acorr_d,1);
    tot_lf8_acorr_d(d,:) = mean(lf8_acorr_d,1);

%     tot_lf5_acorr(d,:) = mean(lf5_acorr,1);
    tot_w8_x(d,:) = mean(w8_x,1);
    tot_w8_x_d(d,:) = mean(w8_x_d,1);
%      tot_w5_x(d,:) = mean(w5_x,1);
%      tot_l85_x(d,:) = mean(l85_x,1);
    
    plot(lags/Fsd,tot_wcv_acorr(d,:),'linewidth',2)
    hold on
    plot(lags/Fsd,tot_lf8_acorr(d,:),'r','linewidth',2)
%     plot(lags/Fsd,tot_lf5_acorr(d,:),'r','linewidth',2)
    legend('MP acorr','LF8 acorr')
    xlim([0 5])
    cell_name = ['L' sess_data(d).layer '_' sess_data(d).cell_type '_' sess_data(d).region '_' sess_data(d).name];
    tname = ['C:\WC_Germany\Cortical_analysis\time_domain\acorrs_' cell_name];
    print('-dpng',tname);
    close

        plot(lags/Fsd,tot_wcv_acorr_d(d,:),'linewidth',2)
    hold on
    plot(lags/Fsd,tot_lf8_acorr_d(d,:),'r','linewidth',2)
%     plot(lags/Fsd,tot_lf5_acorr(d,:),'r','linewidth',2)
    legend('MP acorr','LF8 acorr')
    xlim([0 5])
    cell_name = ['L' sess_data(d).layer '_' sess_data(d).cell_type '_' sess_data(d).region '_' sess_data(d).name];
    tname = ['C:\WC_Germany\Cortical_analysis\time_domain\deriv_acorrs_' cell_name];
    print('-dpng',tname);
    close

    
    plot(lags/Fsd,tot_w8_x(d,:),'linewidth',2)
%     hold on
%         plot(lags/Fsd,tot_w5_x(d,:),'k','linewidth',2)
%     plot(lags/Fsd,tot_l85_x(d,:),'c','linewidth',2)

%     legend('MPxLF8','MPxLF5','LF8xLF5')
    xlim([-1 1]);grid
        cell_name = ['L' sess_data(d).layer '_' sess_data(d).cell_type '_' sess_data(d).region '_' sess_data(d).name];
    tname = ['C:\WC_Germany\Cortical_analysis\time_domain\xcorrs_' cell_name];
    print('-dpng',tname);
    close

        plot(lags/Fsd,tot_w8_x_d(d,:),'linewidth',2)
%     hold on
%         plot(lags/Fsd,tot_w5_x(d,:),'k','linewidth',2)
%     plot(lags/Fsd,tot_l85_x(d,:),'c','linewidth',2)

%     legend('MPxLF8','MPxLF5','LF8xLF5')
    xlim([-0.2 0.2]);grid
        cell_name = ['L' sess_data(d).layer '_' sess_data(d).cell_type '_' sess_data(d).region '_' sess_data(d).name];
    tname = ['C:\WC_Germany\Cortical_analysis\time_domain\deriv_xcorrs_' cell_name];
    print('-dpng',tname);
    close

    
    clear down_w down_8 deriv_w deriv_8 wcv* lf8* lf3* w8_x* w3_x* w3s_x*

end

save C:\WC_Germany\Cortical_analysis\time_domain\time_domain_data tot* Fsd lags