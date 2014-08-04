clear all

load C:\WC_Germany\Persistent_activity\dir_tree_update
load C:\WC_Germany\Persistent_activity\desynch_extract\desynch_points

dsf = 8;
Fsd = 2016/dsf;
minSegLength = 100;
maxLag = 10*Fsd;
niqf = 2016/2;
lcf = .05/niqf;
hcf = 40/niqf;
[b,a] = butter(2,[lcf hcf]);

for d = 1:length(dir_array)

    disp(sprintf('session %d',d))
    cd(dir_array{d});
    pwd

    load used_data wcv_minus_spike lf8

    wcv_f = filtfilt(b,a,wcv_minus_spike);
    lf8_f = filtfilt(b,a,lf8);

    down_w = downsample(wcv_f,dsf);
    down_8 = downsample(lf8,dsf);
    
    %zscore
    down_w = zscore(down_w);
    down_8 = zscore(down_8);

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
            [w8_x(cnt,:),lags] = xcov(down_w,down_8,maxLag,'coeff');
        end
    end

    tot_wcv_acorr(d,:) = mean(wcv_acorr,1);
    tot_lf8_acorr(d,:) = mean(lf8_acorr,1);
    tot_w8_x(d,:) = mean(w8_x,1);

    subplot(2,1,1)
    plot(lags/Fsd,tot_wcv_acorr(d,:),'linewidth',2)
    subplot(2,1,2)
    xlim([0 10])
    plot(lags/Fsd,tot_wcv_acorr(d,:),'linewidth',2)
    xlim([0 2]);grid
    tname = ['C:\WC_Germany\Persistent_activity\time_domain\w_acorr' f_names{d}];
    print('-dpng',tname);
    close

    subplot(2,1,1)
    plot(lags/Fsd,tot_lf8_acorr(d,:),'linewidth',2)
        xlim([0 10])
    subplot(2,1,2)
    plot(lags/Fsd,tot_lf8_acorr(d,:),'linewidth',2)
    xlim([0 2]);grid
    tname = ['C:\WC_Germany\Persistent_activity\time_domain\lf8_acorr' f_names{d}];
    print('-dpng',tname);
    close

    subplot(2,1,1)
    plot(lags/Fsd,tot_w8_x(d,:),'linewidth',2)
    subplot(2,1,2)
    plot(lags/Fsd,tot_w8_x(d,:),'linewidth',2)
    xlim([-2 2]);grid
    tname = ['C:\WC_Germany\Persistent_activity\time_domain\xcor' f_names{d}];
    print('-dpng',tname);
    close

    clear down_w down_8 wcv* lf8

end


save C:\WC_Germany\Persistent_activity\time_domain\td_data tot_wcv_acorr tot_lf8_acorr tot_w8_x