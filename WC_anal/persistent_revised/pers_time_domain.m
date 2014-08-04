clear all
close all

load C:\WC_Germany\persistent_revised\pers_revised_dir
load C:\WC_Germany\persistent_revised\desynch_detect\desynch_times
load C:\WC_Germany\persistent_revised\UDS_dur_raw\UDS_raw_data
load C:\WC_Germany\persistent_revised\UDS_synch_state_dur\UDS_synch_state_dur_data

dsf = 8;
Fsd = 2016/dsf;
minSegLength = 100;
maxLag = 10*Fsd;
niqf = 2016/2;
lcf = .05/niqf;
hcf = 40/niqf;
[b,a] = butter(2,[lcf hcf]);

for d = 1:28

    disp(sprintf('session %d',d))
    cd(dir_array{d});
    pwd

    load used_data wcv_minus_spike lf8 

    wcv_f = filtfilt(b,a,wcv_minus_spike);
    lf8_f = filtfilt(b,a,lf8);
    
    down_w = downsample(wcv_f,dsf);
    down_8 = downsample(lf8_f,dsf);
    
    %zscore
    down_w = zscore(down_w);
    down_8 = zscore(down_8);
    
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
            
            [w8_x(cnt,:),lags] = xcov(down_w,down_8,maxLag,'coeff');
        end
    end

    tot_wcv_acorr(d,:) = mean(wcv_acorr,1);
    tot_lf8_acorr(d,:) = mean(lf8_acorr,1);
    tot_w8_x(d,:) = mean(w8_x,1);
        
    
    plot(lags/Fsd,tot_wcv_acorr(d,:),'linewidth',2)
    hold on
    plot(lags/Fsd,tot_lf8_acorr(d,:),'r','linewidth',2)
    legend('MP acorr','LF8 acorr','LF3 acorr','LF2')
    xlim([0 10])
    tname = ['C:\WC_Germany\persistent_revised\time_domain\acorrs_' f_names{d}];
    print('-dpng',tname);
    close

    plot(lags/Fsd,tot_w8_x(d,:),'linewidth',2)
    legend('MPxLF8')
    xlim([-2 2]);grid
    tname = ['C:\WC_Germany\persistent_revised\time_domain\xcorrs_' f_names{d}];
    print('-dpng',tname);
    close

    clear down_w down_8 wcv* lf8* lf3* w8_x* w3_x* w3s_x*

end

save C:\WC_Germany\persistent_revised\time_domain\time_domain_data tot* Fsd lags