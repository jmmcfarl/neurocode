%% Find persistent activity
clear all
cd C:\WC_Germany\Tetrode_recs\2008-12-16
dsf = 16;
Fsd = 2016/dsf;

windowSize = 60; %in s
windowSlide = 1; %in s
smooth_param = 8; % number of prob bins to smooth
numBins = 200; %number of prob bins

%filter low and medium
niqf = 2016/2;
hcf = 1/niqf;
lcf = .05/niqf;
[b,a] = butter(2,[lcf hcf]);
% [b,a] = butter(2,hcf,'low');
% %set specgram params
% params.fpass = [2 40];
% params.Fs = 2016/4;
% params.tapers = [1 1];
% window = [1 0.1];



    %load data, downsample, and filter
    load lfp_data
    
    lf8 = lf1_data;
    lf15 = lf12_data;

    %     %get specgrams
    %     wcvs = downsample(wcv_minus_spike,4);
    %     lf8s = downsample(lf8,4);
    %     [S8,t8,f8] = mtspecgramc(lf8s,window,params);
    %     [Sw,tw,fw] = mtspecgramc(wcvs,window,params);


    lf8_f = filtfilt(b,a,lf8);
    lf15_f = filtfilt(b,a,lf15);

    lf8_d = downsample(lf8_f,dsf);
    lf15_d = downsample(lf15_f,dsf);

    lf8_d(1:3000) = [];
    lf15_d(1:3000) = [];
    
    lf8_d = zscore(lf8_d);
    lf15_d = zscore(lf15_d);

    total_dur = round(length(lf8_d)/2016*dsf);

    numWins = floor((total_dur-windowSize)/windowSlide);

%     range8 = linspace(min(lf8_d),max(lf8_d),numBins);
%     range15 = linspace(min(lf15_d),max(lf15_d),numBins);
    t_axis = [1:numWins]*windowSlide;
%     diff8 = mean(diff(range8));
%     diff15 = mean(diff(range15));

    %unimodal times
    uni_times_15 = zeros(size(t_axis));
    uni_times_8 = zeros(size(t_axis));
    %     bad_times_w{d} = zeros(size(t_axis{d}));
    %     bad_times_8{d} = zeros(size(t_axis{d}));
    up_peaks_15 = [];
    down_peaks_15 = [];
    up_peaks_8 = [];
    down_peaks_8 = [];
    uni_peaks_15 = zeros(size(t_axis));
    uni_peaks_8 = zeros(size(t_axis));
    uni_heights_15 = zeros(size(t_axis));
    uni_heights_8 = zeros(size(t_axis));
    
    for w = 1:numWins

        begT = (w-1)*windowSlide;
        endT = begT + windowSize;

        begInd = begT*Fsd+1;
        endInd = min((begInd+windowSize*Fsd),length(lf15_d));

        lf8_seg = lf8_d(begInd:endInd);
        lf15_seg = lf15_d(begInd:endInd);

%         lf8_dist(w,:) = jmm_smooth_1d(hist(lf8_seg,range8),smooth_param);
%         lf15_dist(w,:) = jmm_smooth_1d(hist(lf15_seg,range15),smooth_param);
% 
%         lf8_dist(w,:) = lf8_dist(w,:)/sum(lf8_dist(w,:));
%         lf15_dist(w,:) = lf15_dist(w,:)/sum(lf15_dist(w,:));

[lf8_dist(w,:),range8] = gpkde(lf8_seg',-3,[min(lf8_d) max(lf8_d)]);
[lf15_dist(w,:),range15] = gpkde(lf15_seg',-3,[min(lf15_d) max(lf15_d)]);
    diff8 = mean(diff(range8));
    diff15 = mean(diff(range15));


        [peak_heights,peaks_15] = findpeaks(lf15_dist(w,:),'minpeakdistance',round(.5/diff15));

        %get rid of secondary peaks beyond 3 sd above mean
        bad_peaks = find(range15(peaks_15) > 3);
        peaks_15(bad_peaks) = [];
        peak_heights(bad_peaks) = [];

        %if the local distribution is bimodal
        if length(peaks_15) > 1
            
            [dummy,peak_order] = sort(peak_heights,'descend');
            peaks_15 = peaks_15(peak_order);
            up15 = max([peaks_15(1) peaks_15(2)]);
            down15 = min([peaks_15(1) peaks_15(2)]);
            up_peaks_15 = [up_peaks_15 up15];
            down_peaks_15 = [down_peaks_15 down15];
            threshold_15(w) = range15(round(down15*.5 + up15*.5));
        
        
        else  %if the distribution is unimodal
                
            uni_times_15(w) = 1;
               %find skewness of the distribution
               skew_15 = skewness(lf15_seg);
               
               %if the distribution is positively skewed
               if skew_15 > 0
                   
                   %set threshold as first inflection point after the peak
                   [ddif_pval,first_ddif_peak] = findpeaks(diff(lf15_dist(w,peaks_15:end),2),'npeaks',1);
                   threshold_15(w) = range15(peaks_15-1+first_ddif_peak);
                  
                   %if the distribution is negatively skewed
               else
                   
                   [ddif_pval,first_ddif_peak] = findpeaks(diff(lf15_dist(w,1:peaks_15-3),2));
                   threshold_15(w) = range15(first_ddif_peak(end));
                   
               end

        end

        
        [peak_heights,peaks_8] = findpeaks(lf8_dist(w,:),'minpeakdistance',round(.5/diff8));

        bad_peaks = find(range8(peaks_8) > 3);
        peaks_8(bad_peaks) = [];
        peak_heights(bad_peaks) = [];

        if length(peaks_8) > 1
            [dummy,peak_order] = sort(peak_heights,'descend');
            peaks_8 = peaks_8(peak_order);
            up8 = max([peaks_8(1) peaks_8(2)]);
            down8 = min([peaks_8(1) peaks_8(2)]);
            up_peaks_8 = [up_peaks_8 up8];
            down_peaks_8 = [down_peaks_8 down8];
            threshold_8(w) = range8(round(down8*0.5 + up8*0.5));
        else

            uni_times_8(w) = 1;
               %find skewness of the distribution
               skew_8 = skewness(lf8_seg);
               
               %if the distribution is positively skewed
               if skew_8 > 0
                   
                   %set threshold as first inflection point after the peak
                   [ddif_pval,first_ddif_peak] = findpeaks(diff(lf8_dist(w,peaks_8:end),2),'npeaks',1);
                   threshold_8(w) = range8(peaks_8-1+first_ddif_peak);
                  
                   %if the distribution is negatively skewed
               else
                   
                   [ddif_pval,first_ddif_peak] = findpeaks(diff(lf8_dist(w,1:peaks_8),2));
                   threshold_8(w) = range8(first_ddif_peak(end));
                   
               end



        end



    end
    
    



    sm_threshold_15 = jmm_smooth_1d(threshold_15,10);
    sm_threshold_8 = jmm_smooth_1d(threshold_8,10);
  
    Fig = figure(1)
    clf
    set(Fig,'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [40 40]); % paper size is in [width height] format
    set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
    subplot(2,1,1)
    pcolor(t_axis,range15,lf15_dist');shading flat
    hold on
    plot(t_axis,threshold_15,'w.','linewidth',3)
    plot(t_axis,sm_threshold_15,'w','linewidth',3)
    plot(t_axis(find(uni_times_15)),threshold_15(find(uni_times_15)),'r.','linewidth',3)
    plot(t_axis(find(uni_times_15)),sm_threshold_15(find(uni_times_15)),'r.','linewidth',3)
    caxis([0 0.5])

    subplot(2,1,2)
    pcolor(t_axis,range8,lf8_dist');shading flat
    hold on
    plot(t_axis,threshold_8,'w.','linewidth',3)
    plot(t_axis,sm_threshold_8,'w','linewidth',3)
    plot(t_axis(find(uni_times_8)),threshold_8(find(uni_times_8)),'r.','linewidth',3)
    plot(t_axis(find(uni_times_8)),sm_threshold_8(find(uni_times_8)),'r.','linewidth',3)
    caxis([0 0.5])

    t_name = 'C:\WC_Germany\Tetrode_recs\2008-12-16';
    print('-dpng',t_name)
    close all

save up_per_data threshold* 

%% PLOT UP and DOwn State Durations
% %%session 1
% clear all
% load C:\WC_Germany\JMM_Analysis_pyr\sim_record\up_per_data1
% 
%     [up_dist,up_r] = log_hist(up_state_dur,[0.3 50],100);
% [down_dist,down_r] = log_hist(down_state_dur,[0.3 50],100);
% [up_dist8,up_r] = log_hist(up_state_dur8,[0.3 50],100);
% [down_dist8,down_r]  = log_hist(down_state_dur8,[0.3 50],100);
% 


