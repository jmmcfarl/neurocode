clear all

load E:\WC_Germany\JMM_Analysis\file_directories_jmm.mat

window = 30; %window size in sec
win_slide = 2; %window slide

dsf = 8;
Fsd = 2016/dsf;
niqf = Fsd/2;
hif = 2/niqf;
[b,a] = butter(2,hif,'low');

for d = 26:length(dir_array)

    cd(dir_array{d})
    pwd

    load used_data lf8 wcv_minus_spike

    lf8_down = downsample(lf8,dsf);
    wcv_down = downsample(wcv_minus_spike,dsf);

    lf8_down = filtfilt(b,a,lf8_down);
    wcv_down = filtfilt(b,a,wcv_down);

    lf8_down_z = lf8_down - mean(lf8_down);
    lf8_down_z = lf8_down_z/std(lf8_down_z);

    wcv_down_z = wcv_down - mean(wcv_down);
    wcv_down_z = wcv_down_z/std(wcv_down_z);

    record_dur = length(lf8_down_z)/Fsd;
    numWins = floor((record_dur-window)/win_slide);

    for i = 1:numWins

        begPt = (i-1)*win_slide*Fsd+1;
        endPt = begPt+window*Fsd;

        corTime(i) = (begPt+endPt)/2/Fsd;

        lf8Seg = lf8_down_z(begPt:endPt);
        wcvSeg = wcv_down_z(begPt:endPt);

        diff8Seg = [0;diff(lf8Seg)];
        diff8Seg = fgsmooth(diff8Seg,3);
        diffWSeg = [0;diff(wcvSeg)];
        diffWSeg = fgsmooth(diffWSeg,3);

        clf
        subplot(2,1,1)
        plot(lf8Seg)
        hold on
        plot(wcvSeg,'k')
        subplot(2,1,2)
        plot(diff8Seg)
        hold on
        plot(diffWSeg,'k')


    end



end

