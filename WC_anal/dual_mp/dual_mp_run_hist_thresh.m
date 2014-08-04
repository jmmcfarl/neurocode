%% Find persistent activity
clear all
cd C:\WC_Germany\DUAL-MP_recordings
load dual_mp_dir

dsf = 8;
Fsd = 2016/dsf;

windowSize = 60; %in s
windowSlide = 1; %in s
% smooth_param = 15; % number of prob bins to smooth
% numBins = 200; %number of prob bins

%filter low and medium
niqf = 2016/2;
hcf = 1/niqf;
lcf = .05/niqf;
[b,a] = butter(2,[lcf hcf]);

for d = 1:length(dir_array)

    cd(dir_array{d})
    pwd

    %load data, downsample, and filter
    load used_data lf*

    lf8_f = filtfilt(b,a,lf8);
    lf8_d = downsample(lf8_f,dsf);
    lf8_d = zscore(lf8_d);

    lf13_f = filtfilt(b,a,lf13);
    lf13_d = downsample(lf13_f,dsf);
    lf13_d = zscore(lf13_d);

    lf15_f = filtfilt(b,a,lf15);
    lf15_d = downsample(lf15_f,dsf);
    lf15_d = zscore(lf15_d);

    lf16_f = filtfilt(b,a,lf16);
    lf16_d = downsample(lf16_f,dsf);
    lf16_d = zscore(lf16_d);

    if exist('lf14')
        lf14_f = filtfilt(b,a,lf14);
        lf14_d = downsample(lf14_f,dsf);
        lf14_d = zscore(lf14_d);
    end

    %load time data
    total_dur = length(lf8)/2016;
    numWins = floor((total_dur-windowSize)/windowSlide);

    t_axis{d} = (1:numWins)*windowSlide;

    %unimodal times
    uni_times_8{d} = zeros(size(t_axis{d}));
    uni_times_13{d} = zeros(size(t_axis{d}));
    uni_times_15{d} = zeros(size(t_axis{d}));
    uni_times_16{d} = zeros(size(t_axis{d}));
    if exist('lf14')
        uni_times_14{d} = zeros(size(t_axis{d}));
    end

    %initialize peak vectors
    up_peaks_8 = [];
    down_peaks_8 = [];
    uni_peaks_8 = zeros(size(t_axis{d}));
    uni_heights_8 = zeros(size(t_axis{d}));
    up_peaks_13 = [];
    down_peaks_13 = [];
    uni_peaks_13 = zeros(size(t_axis{d}));
    uni_heights_13 = zeros(size(t_axis{d}));
    up_peaks_15 = [];
    down_peaks_15 = [];
    uni_peaks_15 = zeros(size(t_axis{d}));
    uni_heights_15 = zeros(size(t_axis{d}));
    up_peaks_16 = [];
    down_peaks_16 = [];
    uni_peaks_16 = zeros(size(t_axis{d}));
    uni_heights_16 = zeros(size(t_axis{d}));
    if exist('lf14')
        up_peaks_14 = [];
        down_peaks_14 = [];
        uni_peaks_14 = zeros(size(t_axis{d}));
        uni_heights_14 = zeros(size(t_axis{d}));
    end

    for w = 1:numWins

        begT = (w-1)*windowSlide;
        endT = begT + windowSize;

        begInd = begT*Fsd+1;
        endInd = min((begInd+windowSize*Fsd),length(lf8_d));

        lf8_seg = lf8_d(begInd:endInd);
        lf13_seg = lf13_d(begInd:endInd);
        lf15_seg = lf15_d(begInd:endInd);
        lf16_seg = lf16_d(begInd:endInd);
        if exist('lf14')
            lf14_seg = lf14_d(begInd:endInd);
        end

        [lf8_dist{d}(w,:),range8{d}] = gpkde(lf8_seg,-3,[min(lf8_d);max(lf8_d)]);
        [lf13_dist{d}(w,:),range13{d}] = gpkde(lf13_seg,-3,[min(lf13_d);max(lf13_d)]);
        [lf15_dist{d}(w,:),range15{d}] = gpkde(lf15_seg,-3,[min(lf15_d);max(lf15_d)]);
        [lf16_dist{d}(w,:),range16{d}] = gpkde(lf16_seg,-3,[min(lf16_d);max(lf16_d)]);
        if exist('lf14')
            [lf14_dist{d}(w,:),range14{d}] = gpkde(lf14_seg,-3,[min(lf14_d);max(lf14_d)]);
        end

        diff8 = mean(diff(range8{d}));
        diff13 = mean(diff(range13{d}));
        diff15 = mean(diff(range15{d}));
        diff16 = mean(diff(range16{d}));
        if exist('lf14')
            diff14 = mean(diff(range14{d}));
        end

        %% LF8
        [peak_heights,peaks_8] = findpeaks(lf8_dist{d}(w,:),'minpeakdistance',round(.5/diff8));

        %         bad_peaks = find(range8{d}(peaks_8) > 3);
        %         peaks_8(bad_peaks) = [];
        %         peak_heights(bad_peaks) = [];

        if length(peaks_8) > 1
            [dummy,peak_order] = sort(peak_heights,'descend');
            peaks_8 = peaks_8(peak_order);
            up8 = max([peaks_8(1) peaks_8(2)]);
            down8 = min([peaks_8(1) peaks_8(2)]);
            up_peaks_8 = [up_peaks_8 up8];
            down_peaks_8 = [down_peaks_8 down8];
            threshold_8{d}(w) = range8{d}(round(down8*0.5 + up8*0.5));
        else

            uni_times_8{d}(w) = 1;
            %find skewness of the distribution
            skew_8 = skewness(lf8_seg);

            %if the distribution is positively skewed
            if skew_8 > 0

                %set threshold as first inflection point after the peak
                [ddif_pval,first_ddif_peak] = findpeaks(diff(lf8_dist{d}(w,peaks_8:end),2),'npeaks',1);
                threshold_8{d}(w) = range8{d}(peaks_8-1+first_ddif_peak);

                %if the distribution is negatively skewed
            else

                [ddif_pval,first_ddif_peak] = findpeaks(diff(lf8_dist{d}(w,1:peaks_8),2));
                threshold_8{d}(w) = range8{d}(first_ddif_peak(end));

            end

        end

        %% LF13
        [peak_heights,peaks_13] = findpeaks(lf13_dist{d}(w,:),'minpeakdistance',round(.5/diff13));

        if length(peaks_13) > 1
            [dummy,peak_order] = sort(peak_heights,'descend');
            peaks_13 = peaks_13(peak_order);
            up13 = max([peaks_13(1) peaks_13(2)]);
            down13 = min([peaks_13(1) peaks_13(2)]);
            up_peaks_13 = [up_peaks_13 up13];
            down_peaks_13 = [down_peaks_13 down13];
            threshold_13{d}(w) = range13{d}(round(down13*0.5 + up13*0.5));
        else

            uni_times_13{d}(w) = 1;
            %find skewness of the distribution
            skew_13 = skewness(lf13_seg);

            %if the distribution is positively skewed
            if skew_13 > 0

                %set threshold as first inflection point after the peak
                [ddif_pval,first_ddif_peak] = findpeaks(diff(lf13_dist{d}(w,peaks_13:end),2),'npeaks',1);
                threshold_13{d}(w) = range13{d}(peaks_13-1+first_ddif_peak);

                %if the distribution is negatively skewed
            else

                [ddif_pval,first_ddif_peak] = findpeaks(diff(lf13_dist{d}(w,1:peaks_13),2));
                threshold_13{d}(w) = range13{d}(first_ddif_peak(end));

            end

        end

        %% LF15
        [peak_heights,peaks_15] = findpeaks(lf15_dist{d}(w,:),'minpeakdistance',round(.5/diff15));

        %         bad_peaks = find(range8{d}(peaks_8) > 3);
        %         peaks_8(bad_peaks) = [];
        %         peak_heights(bad_peaks) = [];

        if length(peaks_15) > 1
            [dummy,peak_order] = sort(peak_heights,'descend');
            peaks_15 = peaks_15(peak_order);
            up15 = max([peaks_15(1) peaks_15(2)]);
            down15 = min([peaks_15(1) peaks_15(2)]);
            up_peaks_15 = [up_peaks_15 up15];
            down_peaks_15 = [down_peaks_15 down15];
            threshold_15{d}(w) = range15{d}(round(down15*0.5 + up15*0.5));
        else

            uni_times_15{d}(w) = 1;
            %find skewness of the distribution
            skew_15 = skewness(lf15_seg);

            %if the distribution is positively skewed
            if skew_15 > 0

                %set threshold as first inflection point after the peak
                [ddif_pval,first_ddif_peak] = findpeaks(diff(lf15_dist{d}(w,peaks_15:end),2),'npeaks',1);
                threshold_15{d}(w) = range15{d}(peaks_15-1+first_ddif_peak);

                %if the distribution is negatively skewed
            else

                [ddif_pval,first_ddif_peak] = findpeaks(diff(lf15_dist{d}(w,1:peaks_15),2));
                threshold_15{d}(w) = range15{d}(first_ddif_peak(end));

            end

        end
        %% LF16

        [peak_heights,peaks_16] = findpeaks(lf16_dist{d}(w,:),'minpeakdistance',round(.5/diff16));


        if length(peaks_16) > 1
            [dummy,peak_order] = sort(peak_heights,'descend');
            peaks_16 = peaks_16(peak_order);
            up16 = max([peaks_16(1) peaks_16(2)]);
            down16 = min([peaks_16(1) peaks_16(2)]);
            up_peaks_16 = [up_peaks_16 up16];
            down_peaks_16 = [down_peaks_16 down16];
            threshold_16{d}(w) = range16{d}(round(down16*0.5 + up16*0.5));
        else

            uni_times_16{d}(w) = 1;
            %find skewness of the distribution
            skew_16 = skewness(lf16_seg);

            %if the distribution is positively skewed
            if skew_16 > 0

                %set threshold as first inflection point after the peak
                [ddif_pval,first_ddif_peak] = findpeaks(diff(lf16_dist{d}(w,peaks_16:end),2),'npeaks',1);
                threshold_16{d}(w) = range16{d}(peaks_16-1+first_ddif_peak);

                %if the distribution is negatively skewed
            else

                [ddif_pval,first_ddif_peak] = findpeaks(diff(lf16_dist{d}(w,1:peaks_16),2));
                threshold_16{d}(w) = range16{d}(first_ddif_peak(end));

            end

        end

        %%
        if exist('lf14')
            [peak_heights,peaks_14] = findpeaks(lf14_dist{d}(w,:),'minpeakdistance',round(.5/diff14));

            %         bad_peaks = find(range14{d}(peaks_14) > 3);
            %         peaks_14(bad_peaks) = [];
            %         peak_heights(bad_peaks) = [];

            if length(peaks_14) > 1
                [dummy,peak_order] = sort(peak_heights,'descend');
                peaks_14 = peaks_14(peak_order);
                up14 = max([peaks_14(1) peaks_14(2)]);
                down14 = min([peaks_14(1) peaks_14(2)]);
                up_peaks_14 = [up_peaks_14 up14];
                down_peaks_14 = [down_peaks_14 down14];
                threshold_14{d}(w) = range14{d}(round(down14*0.5 + up14*0.5));
            else

                uni_times_14{d}(w) = 1;
                %find skewness of the distribution
                skew_14 = skewness(lf14_seg);

                %if the distribution is positively skewed
                if skew_14 > 0

                    %set threshold as first inflection point after the peak
                    [ddif_pval,first_ddif_peak] = findpeaks(diff(lf14_dist{d}(w,peaks_14:end),2),'npeaks',1);
                    threshold_14{d}(w) = range14{d}(peaks_14-1+first_ddif_peak);

                    %if the distribution is negatively skewed
                else

                    [ddif_pval,first_ddif_peak] = findpeaks(diff(lf14_dist{d}(w,1:peaks_14),2));
                    threshold_14{d}(w) = range14{d}(first_ddif_peak(end));

                end

            end
        end

        %%


    end

    
    niqf = 1/2;
    [b2,a2] = butter(2,.01/niqf,'low');
    %     sm_threshold_w{d} = jmm_smooth_1d_cor(threshold_w{d},10);
    %     sm_threshold_8{d} = jmm_smooth_1d_cor(threshold_8{d},10);

    sm_threshold_8{d} = filtfilt(b2,a2,threshold_8{d});
    sm_threshold_13{d} = filtfilt(b2,a2,threshold_13{d});
    sm_threshold_15{d} = filtfilt(b2,a2,threshold_15{d});
    sm_threshold_16{d} = filtfilt(b2,a2,threshold_16{d});
if exist('lf14')
       sm_threshold_14{d} = filtfilt(b2,a2,threshold_14{d});
end

    Fig = figure(1)
    clf
    set(Fig,'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [40 40]);% paper size is in [width height] format
    set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
    subplot(2,1,1)
    imagesc(t_axis{d},range8{d},lf8_dist{d}');shading flat
    hold on
    plot(t_axis{d},sm_threshold_8{d},'w','linewidth',1)
    plot(t_axis{d}(find(uni_times_8{d})),sm_threshold_8{d}(find(uni_times_8{d})),'r.','MarkerSize',15)
    caxis([0 0.7]);colorbar
    subplot(2,1,2)
    imagesc(t_axis{d},range13{d},lf13_dist{d}');shading flat
    hold on
    plot(t_axis{d},sm_threshold_13{d},'w','linewidth',1)
    plot(t_axis{d}(find(uni_times_13{d})),sm_threshold_13{d}(find(uni_times_13{d})),'r.','MarkerSize',15)
    caxis([0 0.7]);colorbar
    t_name = ['C:\WC_Germany\DUAL-MP_recordings\run_hist_thresh\lf13_' f_names{d}];
    print('-dpng',t_name)
    close all

    Fig = figure(1)
    clf
    set(Fig,'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [40 40]);% paper size is in [width height] format
    set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
    subplot(2,1,1)
    imagesc(t_axis{d},range8{d},lf8_dist{d}');shading flat
    hold on
    plot(t_axis{d},sm_threshold_8{d},'w','linewidth',1)
    plot(t_axis{d}(find(uni_times_8{d})),sm_threshold_8{d}(find(uni_times_8{d})),'r.','MarkerSize',15)
    caxis([0 0.7]);colorbar
    subplot(2,1,2)
    imagesc(t_axis{d},range15{d},lf15_dist{d}');shading flat
    hold on
    plot(t_axis{d},sm_threshold_15{d},'w','linewidth',1)
    plot(t_axis{d}(find(uni_times_15{d})),sm_threshold_15{d}(find(uni_times_15{d})),'r.','MarkerSize',15)
    caxis([0 0.7]);colorbar
    t_name = ['C:\WC_Germany\DUAL-MP_recordings\run_hist_thresh\lf15_' f_names{d}];
    print('-dpng',t_name)
    close all

        Fig = figure(1)
    clf
    set(Fig,'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [40 40]);% paper size is in [width height] format
    set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
    subplot(2,1,1)
    imagesc(t_axis{d},range8{d},lf8_dist{d}');shading flat
    hold on
    plot(t_axis{d},sm_threshold_8{d},'w','linewidth',1)
    plot(t_axis{d}(find(uni_times_8{d})),sm_threshold_8{d}(find(uni_times_8{d})),'r.','MarkerSize',15)
    caxis([0 0.7]);colorbar
    subplot(2,1,2)
    imagesc(t_axis{d},range16{d},lf16_dist{d}');shading flat
    hold on
    plot(t_axis{d},sm_threshold_16{d},'w','linewidth',1)
    plot(t_axis{d}(find(uni_times_16{d})),sm_threshold_16{d}(find(uni_times_16{d})),'r.','MarkerSize',15)
    caxis([0 0.7]);colorbar
    t_name = ['C:\WC_Germany\DUAL-MP_recordings\run_hist_thresh\lf16_' f_names{d}];
    print('-dpng',t_name)
    close all

    if exist('lf14')
           Fig = figure(1)
    clf
    set(Fig,'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [40 40]);% paper size is in [width height] format
    set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
    subplot(2,1,1)
    imagesc(t_axis{d},range8{d},lf8_dist{d}');shading flat
    hold on
    plot(t_axis{d},sm_threshold_8{d},'w','linewidth',1)
    plot(t_axis{d}(find(uni_times_8{d})),sm_threshold_8{d}(find(uni_times_8{d})),'r.','MarkerSize',15)
    caxis([0 0.7]);colorbar
    subplot(2,1,2)
    imagesc(t_axis{d},range14{d},lf14_dist{d}');shading flat
    hold on
    plot(t_axis{d},sm_threshold_14{d},'w','linewidth',1)
    plot(t_axis{d}(find(uni_times_14{d})),sm_threshold_14{d}(find(uni_times_14{d})),'r.','MarkerSize',15)
    caxis([0 0.7]);colorbar
    t_name = ['C:\WC_Germany\DUAL-MP_recordings\run_hist_thresh\lf14_' f_names{d}];
    print('-dpng',t_name)
    close all
 
    end
end

save C:\WC_Germany\DUAL-MP_recordings\run_hist_thresh\UDS_dur_data threshold*
