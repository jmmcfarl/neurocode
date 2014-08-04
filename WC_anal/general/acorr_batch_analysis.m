%% GET UP STATE AMPLITUDE BATCH

clear all

load E:\WC_Germany\JMM_Analysis\file_directories_jmm.mat


%initialize overall variables

d = 1;

while d <= length(dir_array)

    cd(dir_array{d})
    pwd

    load used_data wcv_minus_spike lf8 CSC8_SampleFrequencies
    load used_data_multLFP lf2 lf5
    load used_data_multLFP_2 lf3 lf7

    Fs = mean(CSC8_SampleFrequencies);
    maxLag = 10*Fs; %max lag in sec


    
    %filter data
    nyqf = Fs/2;
    lif = 0.05/nyqf;
    hif = 2/nyqf;
    [b,a] = butter(2, [lif hif]);
    lofwcv = filtfilt(b,a,wcv_minus_spike);
    loflf2 = filtfilt(b,a,lf2);
    loflf3 = filtfilt(b,a,lf3);
    loflf5 = filtfilt(b,a,lf5);
%     loflf7 = filtfilt(b,a,lf7);
    loflf8 = filtfilt(b,a,lf8);


    windowSize = 50; %window size in s
    windowSizeBins = round(windowSize*Fs);
    timeSlide = 25; %amount to slide time window in s
    binSlide = round(timeSlide*Fs);
    numWins = round((length(lf3)-windowSizeBins)/binSlide); %number of windows

    for i = 1:numWins

        windowCent = (i-1)*binSlide+round(windowSizeBins/2);
        windowBeg = windowCent - round(windowSizeBins/2)+1;
        windowEnd = windowCent + round(windowSizeBins/2);

        wcvSeg = lofwcv(windowBeg:windowEnd);
        lf2Seg = loflf2(windowBeg:windowEnd);
        lf3Seg = loflf3(windowBeg:windowEnd);
        lf5Seg = loflf5(windowBeg:windowEnd);
%         lf7Seg = loflf7(windowBeg:windowEnd);
        lf8Seg = loflf8(windowBeg:windowEnd);

        [temp_acorr_LF2(i,:),lagVec] = xcov(lf2Seg,lf2Seg,maxLag,'coef');
        [temp_acorr_LF3(i,:),lagVec] = xcov(lf3Seg,lf3Seg,maxLag,'coef');
        [temp_acorr_LF5(i,:),lagVec] = xcov(lf5Seg,lf5Seg,maxLag,'coef');
%         [temp_acorr_LF7(i,:),lagVec] = xcov(lf7Seg,lf7Seg,maxLag,'coef');
        [temp_acorr_LF8(i,:),lagVec] = xcov(lf8Seg,lf8Seg,maxLag,'coef');
        [temp_acorr_WCV(i,:),lagVec] = xcov(wcvSeg,wcvSeg,maxLag,'coef');
    end
    
    acorr_LF2(d,:) = mean(temp_acorr_LF2);
    acorr_LF3(d,:) = mean(temp_acorr_LF3);
    acorr_LF5(d,:) = mean(temp_acorr_LF5);
%     acorr_LF7(d,:) = mean(temp_acorr_LF7);
    acorr_LF8(d,:) = mean(temp_acorr_LF8);
    acorr_WCV(d,:) = mean(temp_acorr_WCV);


    d = d +1

    save E:\WC_Germany\JMM_Analysis\overall_data_jmm_acorr d acorr* lagVec

    %prepare for next iteration
    clear all

    load E:\WC_Germany\JMM_Analysis\file_directories_jmm.mat
    load E:\WC_Germany\JMM_Analysis\overall_data_jmm_acorr.mat


end
