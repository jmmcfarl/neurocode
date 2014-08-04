%% Find persistent activity
clear all

load C:\WC_Germany\persistent_revised\pers_revised_dir

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

for d = 1:28

    cd(dir_array{d})
    pwd
    
    %load data, downsample, and filter
    load used_data lf8 wcv_minus_spike wcv

    lf8_f = filtfilt(b,a,lf8);
    wcv_f = filtfilt(b,a,wcv_minus_spike);

    lf8_d = downsample(lf8_f,dsf);
    wcv_d = downsample(wcv_f,dsf);

    lf8_d = zscore(lf8_d);
    wcv_d = zscore(wcv_d);

    %load time data
    total_dur = length(wcv)/2016;
    numWins = floor((total_dur-windowSize)/windowSlide);

    t_axis{d} = (1:numWins)*windowSlide;
    
    %unimodal times
    uni_times_w{d} = zeros(size(t_axis{d}));
    uni_times_8{d} = zeros(size(t_axis{d}));

    %initialize peak vectors
    up_peaks_w = [];
    down_peaks_w = [];
    up_peaks_8 = [];
    down_peaks_8 = [];
    uni_peaks_w = zeros(size(t_axis{d}));
    uni_peaks_8 = zeros(size(t_axis{d}));
    uni_heights_w = zeros(size(t_axis{d}));
    uni_heights_8 = zeros(size(t_axis{d}));
    
    for w = 1:numWins

        begT = (w-1)*windowSlide;
        endT = begT + windowSize;

        begInd = begT*Fsd+1;
        endInd = min((begInd+windowSize*Fsd),length(wcv_d));

        lf8_seg = lf8_d(begInd:endInd);
        wcv_seg = wcv_d(begInd:endInd);

        [lf8_dist{d}(w,:),range8{d}] = gpkde(lf8_seg,-3,[min(lf8_d);max(lf8_d)]);
        [wcv_dist{d}(w,:),rangew{d}] = gpkde(wcv_seg,-3,[min(wcv_d);max(wcv_d)]);

        diff8 = mean(diff(range8{d}));
        diffw = mean(diff(rangew{d}));

        %find peaks in the distribution that are separated by at least 0.5
        %SD
        [peak_heights,peaks_w] = findpeaks(wcv_dist{d}(w,:),'minpeakdistance',round(.5/diffw));

        %if the local distribution is bimodal
        if length(peaks_w) > 1
            
            %sort peaks in density according to their mass
            [dummy,peak_order] = sort(peak_heights,'descend');     
            peaks_w = peaks_w(peak_order);
            
            %locate peak corresponding to up state and down state
            upw = max([peaks_w(1) peaks_w(2)]);
            downw = min([peaks_w(1) peaks_w(2)]);
            up_peaks_w = [up_peaks_w upw];
            down_peaks_w = [down_peaks_w downw];
            threshold_w{d}(w) = rangew{d}(round(downw*.5 + upw*.5));
        
        
        else  %if the distribution is unimodal
                
               uni_times_w{d}(w) = 1;
               
               %find skewness of the distribution
               skew_w = skewness(wcv_seg);
               
               %if the distribution is positively skewed
               if skew_w > 0
                   
                   %set threshold as first inflection point after the peak
                   [ddif_pval,first_ddif_peak] = findpeaks(diff(wcv_dist{d}(w,peaks_w:end),2),'npeaks',1);
                   threshold_w{d}(w) = rangew{d}(peaks_w-1+first_ddif_peak);
                  
                   %if the distribution is negatively skewed
               else
                   
                   [ddif_pval,first_ddif_peak] = findpeaks(diff(wcv_dist{d}(w,1:peaks_w-3),2));
                   threshold_w{d}(w) = rangew{d}(first_ddif_peak(end));
                   
               end

        end

        %repeat same process for lfp
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



    end

    niqf = 1/2;
    [b2,a2] = butter(2,.01/niqf,'low');
%     sm_threshold_w{d} = jmm_smooth_1d_cor(threshold_w{d},10);
%     sm_threshold_8{d} = jmm_smooth_1d_cor(threshold_8{d},10);
  
sm_threshold_w{d} = filtfilt(b2,a2,threshold_w{d});
sm_threshold_8{d} = filtfilt(b2,a2,threshold_8{d});

    Fig = figure(1)
    clf
    set(Fig,'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [40 40]);% paper size is in [width height] format
    set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
    subplot(2,1,1)
    imagesc(t_axis{d},rangew{d},wcv_dist{d}');shading flat
    hold on
%     plot(t_axis{d},threshold_w{d},'w.','linewidth',3)
    plot(t_axis{d},sm_threshold_w{d},'w','linewidth',1)
%     plot(t_axis{d}(find(uni_times_w{d})),threshold_w{d}(find(uni_times_w{d})),'r.','linewidth',3)
    plot(t_axis{d}(find(uni_times_w{d})),sm_threshold_w{d}(find(uni_times_w{d})),'r.','MarkerSize',15)
    caxis([0 0.7]);colorbar
    subplot(2,1,2)
    imagesc(t_axis{d},range8{d},lf8_dist{d}');shading flat
    hold on
%     plot(t_axis{d},threshold_8{d},'w.','linewidth',3)
    plot(t_axis{d},sm_threshold_8{d},'w','linewidth',1)
%     plot(t_axis{d}(find(uni_times_8{d})),threshold_8{d}(find(uni_times_8{d})),'r.','linewidth',3)
    plot(t_axis{d}(find(uni_times_8{d})),sm_threshold_8{d}(find(uni_times_8{d})),'r.','MarkerSize',15)
    caxis([0 0.7]);colorbar


    t_name = ['C:\WC_Germany\persistent_revised\run_hist_thresh\' f_names{d}];
    print('-dpng',t_name)
    close all

% save C:\WC_Germany\Persistent_activity\run_hist_thresh\UDS_dur_data_new_method threshold* 

end

save C:\WC_Germany\persistent_revised\run_hist_thresh\UDS_dur_data threshold*
