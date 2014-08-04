clear all
close all
cd C:\WC_Germany\DUAL-MP_recordings
load dual_mp_dir
load desynch_detect\desynch_times
load UDS_dur_raw\UDS_raw_data

dsf = 4;
Fsd = 2016/dsf;
minSegLength = 100;
maxLag = 10*Fsd;
niqf = 2016/2;
lcf = .2/niqf;
hcf = 40/niqf;
[b,a] = butter(2,[lcf hcf]);

for d = 1:length(dir_array)

    disp(sprintf('session %d',d))
    cd(dir_array{d});
    pwd

    load used_data lf*

    lf8_f = filtfilt(b,a,lf8);
    down_8 = downsample(lf8_f,dsf);
    down_8 = zscore(down_8);

        lf13_f = filtfilt(b,a,lf13);
    down_13 = downsample(lf13_f,dsf);
    down_13 = zscore(down_13);

        lf15_f = filtfilt(b,a,lf15);
    down_15 = downsample(lf15_f,dsf);
    down_15 = zscore(down_15);

        lf16_f = filtfilt(b,a,lf16);
    down_16 = downsample(lf16_f,dsf);
    down_16 = zscore(down_16);

    if exist('lf14')
                lf14_f = filtfilt(b,a,lf14);
    down_14 = downsample(lf14_f,dsf);
    down_14 = zscore(down_14);
    end
    
    desynch_start = desynch_start_times{d}*Fsd;
    desynch_stop = desynch_stop_times{d}*Fsd;

    if ~isempty(desynch_start)
        sMarkers = [[1;desynch_stop'] [desynch_start';length(down_8)]]
    else
        sMarkers = [1 length(down_8)];
    end

    cnt = 0;
    for i = 1:size(sMarkers,1)
        if (sMarkers(i,2)-sMarkers(i,1))/Fsd > minSegLength
            cnt = cnt+1;
            [lf8_acorr(cnt,:),lags] = xcov(down_8,maxLag,'coeff');
            [lf13_acorr(cnt,:),lags] = xcov(down_13,maxLag,'coeff');
            [lf15_acorr(cnt,:),lags] = xcov(down_15,maxLag,'coeff');
            [lf16_acorr(cnt,:),lags] = xcov(down_16,maxLag,'coeff');
          if exist('lf14')
                          [lf14_acorr(cnt,:),lags] = xcov(down_14,maxLag,'coeff');
          end
            [x8_13(cnt,:),lags] = xcov(down_8,down_13,maxLag,'coeff');
                        [x8_15(cnt,:),lags] = xcov(down_8,down_15,maxLag,'coeff');
            [x8_16(cnt,:),lags] = xcov(down_8,down_16,maxLag,'coeff');
if exist('lf14')
                [x8_14(cnt,:),lags] = xcov(down_8,down_14,maxLag,'coeff');
end
        end
    end

    tot_lf8_acorr(d,:) = mean(lf8_acorr,1);
        tot_lf13_acorr(d,:) = mean(lf13_acorr,1);
    tot_lf15_acorr(d,:) = mean(lf15_acorr,1);
            tot_lf16_acorr(d,:) = mean(lf16_acorr,1);

    tot_x8_13(d,:) = mean(x8_13,1);
    tot_x8_15(d,:) = mean(x8_15,1);
    tot_x8_16(d,:) = mean(x8_16,1);
       
    if exist('lf14')
                    tot_lf14_acorr(d,:) = mean(lf14_acorr,1);
            tot_x8_14(d,:) = mean(x8_14,1);
    end
    
      cur_m = max(tot_x8_13(d,:));
plot(lags/Fsd,tot_x8_13(d,:),'linewidth',2)
ylim([cur_m - 0.025 cur_m])
        xlim([-0.05 0.05]);grid
    tname = ['C:\WC_Germany\DUAL-MP_recordings\time_domain\xcorrs2_8_13_' f_names{d}];
    print('-dpng',tname), close;

       cur_m = max(tot_x8_15(d,:));
plot(lags/Fsd,tot_x8_15(d,:),'linewidth',2)
ylim([cur_m - 0.025 cur_m])
        xlim([-0.05 0.05]);grid
    tname = ['C:\WC_Germany\DUAL-MP_recordings\time_domain\xcorrs2_8_15_' f_names{d}];
    print('-dpng',tname), close;
   
          cur_m = max(tot_x8_16(d,:));
plot(lags/Fsd,tot_x8_16(d,:),'linewidth',2)
ylim([cur_m - 0.025 cur_m])
        xlim([-0.05 0.05]);grid
    tname = ['C:\WC_Germany\DUAL-MP_recordings\time_domain\xcorrs2_8_16_' f_names{d}];
    print('-dpng',tname), close;

    if exist('lf14')
             cur_m = max(tot_x8_14(d,:));
plot(lags/Fsd,tot_x8_14(d,:),'linewidth',2)
ylim([cur_m - 0.025 cur_m])
        xlim([-0.05 0.05]);grid
    tname = ['C:\WC_Germany\DUAL-MP_recordings\time_domain\xcorrs2_8_14_' f_names{d}];
    print('-dpng',tname), close;
 
    end
    if exist('lf14')
          plot(lags/Fsd,tot_lf8_acorr(d,:),'Color',[0.8 0.8 0.7],'linewidth',2)
    hold on
    plot(lags/Fsd,tot_lf13_acorr(d,:),'linewidth',2)
    plot(lags/Fsd,tot_lf14_acorr(d,:),'c','linewidth',2)
        plot(lags/Fsd,tot_lf15_acorr(d,:),'r','linewidth',2)
        plot(lags/Fsd,tot_lf16_acorr(d,:),'k','linewidth',2)
    legend('LF8','LF13','LF14','LF15','LF16')
    xlim([0 10])
    tname = ['C:\WC_Germany\DUAL-MP_recordings\time_domain\acorrs2_' f_names{d}];
    print('-dpng',tname);
    close

    plot(lags/Fsd,tot_x8_13(d,:),'linewidth',2); hold on; grid
    plot(lags/Fsd,tot_x8_14(d,:),'c','linewidth',2)
        plot(lags/Fsd,tot_x8_15(d,:),'r','linewidth',2)
    plot(lags/Fsd,tot_x8_16(d,:),'k','linewidth',2)

    legend('LF8xLF13','LF8xLF14','LF8xLF15','LF8xLF16')
    xlim([-2 2]);grid
    tname = ['C:\WC_Germany\DUAL-MP_recordings\time_domain\xcorrs2_' f_names{d}];
    print('-dpng',tname);
        xlim([-0.05 0.05]);grid
    tname = ['C:\WC_Germany\DUAL-MP_recordings\time_domain\xcorrs2_z_' f_names{d}];
    print('-dpng',tname);

    close

    
    else
    plot(lags/Fsd,tot_lf8_acorr(d,:),'Color',[0.8 0.8 0.7],'linewidth',2)
    hold on
    plot(lags/Fsd,tot_lf13_acorr(d,:),'linewidth',2)
        plot(lags/Fsd,tot_lf15_acorr(d,:),'r','linewidth',2)
        plot(lags/Fsd,tot_lf16_acorr(d,:),'k','linewidth',2)
    legend('LF8','LF13','LF15','LF16')
    xlim([0 10])
    tname = ['C:\WC_Germany\DUAL-MP_recordings\time_domain\acorrs2_' f_names{d}];
    print('-dpng',tname);
    close

    plot(lags/Fsd,tot_x8_13(d,:),'linewidth',2); hold on; grid
        plot(lags/Fsd,tot_x8_15(d,:),'r','linewidth',2)
    plot(lags/Fsd,tot_x8_16(d,:),'k','linewidth',2)

    legend('LF8xLF13','LF8xLF15','LF8xLF16')
    xlim([-2 2]);grid
    tname = ['C:\WC_Germany\DUAL-MP_recordings\time_domain\xcorrs2_' f_names{d}];
    print('-dpng',tname);
        xlim([-0.05 0.05]);grid
    tname = ['C:\WC_Germany\DUAL-MP_recordings\time_domain\xcorrs2_z_' f_names{d}];
    print('-dpng',tname);

    close
    end
    
    clear down_* lf* x8* 

end

save C:\WC_Germany\DUAL-MP_recordings\time_domain\time_domain_data tot* Fsd lags