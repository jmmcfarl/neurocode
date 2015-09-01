%% GET UP STATE AMPLITUDE BATCH

clear all

load E:\WC_Germany\JMM_Analysis\file_directories_jmm.mat


%initialize overall variables

d = 1;

while d <= length(dir_array)

    cd(dir_array{d})
    pwd

    load used_data wcv_minus_spike lf8 CSC8_SampleFrequencies

    Fs = mean(CSC8_SampleFrequencies);
    

    %filter data
    nyqf = Fs/2;
    hif = 1/nyqf;
    lif = 0.05/nyqf;
    [b,a] = butter(2, [lif hif]);
    lofwcv = filtfilt(b,a,wcv_minus_spike);
    loflf8 = filtfilt(b,a,lf8);

    %zscore data
    lofwcv = lofwcv - mean(lofwcv);
    lofwcv = lofwcv/std(lofwcv);
    loflf8 = loflf8 - mean(loflf8);
    loflf8 = loflf8/std(loflf8);

    windowSize = 50; %window size in s
    windowSizeBins = round(windowSize*Fs);
    timeSlide = 5; %amount to slide time window in s
    binSlide = round(timeSlide*Fs);
    numWins = round((length(lf8)-windowSizeBins)/binSlide); %number of windows

    wcv_thresh = (quantile(lofwcv,0.9)+quantile(lofwcv,0.1))/2;
    lf8_thresh = (quantile(loflf8,0.9)+quantile(loflf8,0.1))/2;
    
    for i = 1:numWins

        windowCent = (i-1)*binSlide+round(windowSizeBins/2);
        windowBeg = windowCent - round(windowSizeBins/2)+1;
        windowEnd = windowCent + round(windowSizeBins/2);

        wcvSeg = lofwcv(windowBeg:windowEnd);
        lf8Seg = loflf8(windowBeg:windowEnd);

        wcv_upTime = length(find(wcvSeg>wcv_thresh));
        wcv_downTime = windowSizeBins - wcv_upTime;
        wcv_DCI{d}(i) = (wcv_upTime - wcv_downTime)/(wcv_upTime + wcv_downTime);
        
        lf8_upTime = length(find(lf8Seg>lf8_thresh));
        lf8_downTime = windowSizeBins - lf8_upTime;
        lf8_DCI{d}(i) = (lf8_upTime - lf8_downTime)/(lf8_upTime + lf8_downTime);
    
    end
    

    std_lf8_DCI = std(lf8_DCI{d});
    upperCI_lf8 = lf8_DCI{d}+2*std_lf8_DCI;
    lowerCI_lf8 = lf8_DCI{d}-2*std_lf8_DCI;
    
   plot(wcv_DCI{d},'k')
   hold on
   plot(upperCI_lf8,'r')
   plot(lowerCI_lf8,'r')
   plot(lf8_DCI{d})
   tname = ['E:\WC_Germany\JMM_Analysis\newDCI\CI\' f_names{d}];
   print('-dpng',tname);
   close all
   
    d = d +1

    save E:\WC_Germany\JMM_Analysis\overall_data_jmm_newDCI_test d *DCI

    %prepare for next iteration
    clear all

    load E:\WC_Germany\JMM_Analysis\file_directories_jmm.mat
    load E:\WC_Germany\JMM_Analysis\overall_data_jmm_newDCI_test.mat


end
