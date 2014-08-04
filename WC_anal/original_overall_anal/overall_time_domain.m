clear all
close all

load('C:\WC_Germany\overall_calcs\overall_dir.mat')
load C:\WC_Germany\overall_calcs\desynch_extract\desynch_points

dsf = 8;
Fsd = 2016/dsf;
minSegLength = 100;
maxLag = 10*Fsd;
niqf = 2016/2;
lcf = .05/niqf;
hcf = 40/niqf;
[b,a] = butter(2,[lcf hcf]);

for d = 1:length(over_dir)

    disp(sprintf('session %d',d))
    cd(over_dir{d});
    pwd

    load used_data wcv_minus_spike lf8 lf3 lf2

    wcv_f = filtfilt(b,a,wcv_minus_spike);
    lf2_f = filtfilt(b,a,lf2);
    lf8_f = filtfilt(b,a,lf8);
    lf3_f = filtfilt(b,a,lf3);
    
    down_w = downsample(wcv_f,dsf);
    down_2 = downsample(lf2_f,dsf);
    down_8 = downsample(lf8_f,dsf);
    down_3 = downsample(lf3_f,dsf);
    
    %zscore
    down_w = zscore(down_w);
    down_8 = zscore(down_8);
    down_3 = zscore(down_3);
    
%     % minimize cortical correlation with hipp lfp
%     temp_cor = corrcoef(down_3,down_8);
%     hc_cor = temp_cor(2,1)
%     down_3s = down_3-hc_cor*down_8;

    if ~isempty(desynch_start{d})
        sMarkers = [[1;desynch_stop{d}'] [desynch_start{d}';length(down_w)]]
    else
        sMarkers = [1 length(down_w)];
    end

    cnt = 0;
    for i = 1:size(sMarkers,1)
        if (sMarkers(i,2)-sMarkers(i,1))/Fsd > minSegLength
            cnt = cnt+1;
            
            [wcv_acorr(cnt,:),lags] = xcov(down_w,maxLag,'coeff');
            [lf8_acorr(cnt,:),lags] = xcov(down_8,maxLag,'coeff');
            [lf3_acorr(cnt,:),lags] = xcov(down_3,maxLag,'coeff');
%             [lf3s_acorr(cnt,:),lags] = xcov(down_3s,maxLag,'coeff');
            [lf2_acorr(cnt,:),lags] = xcov(down_2,maxLag,'coeff');
            
            [w8_x(cnt,:),lags] = xcov(down_w,down_8,maxLag,'coeff');
            [w3_x(cnt,:),lags] = xcov(down_w,down_3,maxLag,'coeff');
%             [w3s_x(cnt,:),lags] = xcov(down_w,down_3s,maxLag,'coeff');
%             [x3_8(cnt,:),lags] = xcov(down_3,down_8,maxlag,'coeff');
            [w2_x(cnt,:),lags] = xcov(down_w,down_2,maxLag,'coeff');
        end
    end

    tot_wcv_acorr(d,:) = mean(wcv_acorr,1);
    tot_lf8_acorr(d,:) = mean(lf8_acorr,1);
    tot_w8_x(d,:) = mean(w8_x,1);
    tot_w2_x(d,:) = mean(w2_x,1);
    
    tot_lf3_acorr(d,:) = mean(lf3_acorr,1);
    tot_lf2_acorr(d,:) = mean(lf2_acorr,1);
%     tot_lf3s_acorr(d,:) = mean(lf3s_acorr,1);
    tot_w3_x(d,:) = mean(w3_x,1);
%     tot_w3s_x(d,:) = mean(w3s_x,1);
    
    
    plot(lags/Fsd,tot_wcv_acorr(d,:),'linewidth',2)
    hold on
    plot(lags/Fsd,tot_lf8_acorr(d,:),'r','linewidth',2)
    plot(lags/Fsd,tot_lf3_acorr(d,:),'g','linewidth',2)
    plot(lags/Fsd,tot_lf2_acorr(d,:),'k','linewidth',2)
    legend('MP acorr','LF8 acorr','LF3 acorr','LF2')
    xlim([0 10])
    tname = ['C:\WC_Germany\overall_calcs\time_domain\all_acorrs' num2str(cell_type(d)) '_' over_names{d}];
    print('-dpng',tname);
    close

    plot(lags/Fsd,tot_w8_x(d,:),'linewidth',2)
    hold on
    plot(lags/Fsd,tot_w3_x(d,:),'g','linewidth',2)
    plot(lags/Fsd,tot_w2_x(d,:),'k','linewidth',2)
%     plot(lags/Fsd,tot_w3s_x(d,:),'k','linewidth',2)
    legend('MPxLF8','MPxLF3','MPxLF2')
    xlim([-2 2]);grid
    tname = ['C:\WC_Germany\overall_calcs\time_domain\all_xcorrs' num2str(cell_type(d)) '_' over_names{d}];
    print('-dpng',tname);
    close

    clear down_w down_8 wcv* lf8* lf3* w8_x* w3_x* w3s_x*

end

save C:\WC_Germany\overall_calcs\time_domain\time_domain_data tot* Fsd lags