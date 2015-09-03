%% UDS transition Detector V2

clear all

dir_array{1} = 'E:\WC_Germany\Entorhinal-WC_LFP\2007-05-24_CWC_LFP_B\2007-5-24_16-21-57\part1';
dir_array{2} = 'E:\WC_Germany\Entorhinal-WC_LFP\2007-05-24_CWC_LFP_B\2007-5-24_16-21-57\part2';
dir_array{3} = 'E:\WC_Germany\Entorhinal-WC_LFP\2007-05-25_CWC_LFP_B\2007-5-25_20-24-13\part1';
dir_array{4} = 'E:\WC_Germany\Entorhinal-WC_LFP\2007-05-25_CWC_LFP_B\2007-5-25_20-24-13\part2';
dir_array{5} = 'E:\WC_Germany\Entorhinal-WC_LFP\2007-05-31_CWC_LFP\2007-5-31_17-59-31';
dir_array{6} = 'E:\WC_Germany\Entorhinal-WC_LFP\2007-06-03_CWC_LFP_B\2007-6-3_20-29-2';
dir_array{7} = 'E:\WC_Germany\Entorhinal-WC_LFP\2007-06-04_CWC_LFP_B\2007-6-4_16-24-12';
dir_array{8} = 'E:\WC_Germany\Entorhinal-WC_LFP\2006-04-05_CWC_LFP_A\2006-4-5_19-10-27';
dir_array{9} = 'E:\WC_Germany\Entorhinal-WC_LFP\2007-05-24_CWC_LFP_A\2007-5-24_12-48-28';
dir_array{10} = 'E:\WC_Germany\Entorhinal-WC_LFP\2007-05-28_CWC_LFP_A\2007-5-28_18-27-1';
dir_array{11} = 'E:\WC_Germany\Entorhinal-WC_LFP\2007-06-03_CWC_LFP_A\2007-6-3_19-1-5';
dir_array{12} = 'E:\WC_Germany\Entorhinal-WC_LFP\2007-06-26_CWC_LFP_B\2007-06-26_CWC_LFP_B\2007-6-26_18-47-13\part1';
dir_array{13} = 'E:\WC_Germany\Entorhinal-WC_LFP\2007-06-26_CWC_LFP_B\2007-06-26_CWC_LFP_B\2007-6-26_18-47-13\part2';

f_names{1} = 'pyr_cl_2007_05_24_B_p1';
f_names{2} = 'pyr_cl_2007_05_24_B_p2';
f_names{3} = 'pyr_cl_2007_05_25_B_p1';
f_names{4} = 'pyr_cl_2007_05_25_B_p2';
f_names{5} = 'pyr_cl_2007_05_31';
f_names{6} = 'pyr_cl_2007_06_03_B';
f_names{7} = 'pyr_cl_2007_06_04_B';
f_names{8} = 'stell_cl_2006_04_05_A';
f_names{9} = 'stell_cl_2007_05_24_A';
f_names{10} = 'stell_cl_2007_05_28_A';
f_names{11} = 'stell_cl_2007_06_03_A';
f_names{12} = 'stell_cl_2007_06_26_B_p1';
f_names{13} = 'stell_cl_2007_06_26_B_p2';

lf8_period_overall = [];
wcv_period_overall = [];

for d = 12

    cd(dir_array{d})
    pwd

    if exist('all_eeg_data.mat')
        load all_eeg_data
    elseif exist('part1_eeg_data.mat')
        load part1_eeg_data
    else
        load part2_eeg_data
    end

    clear CSC2* CSC3* CSC4* CSC5* CSC6* CSC7*
    clear *NumberValidSamples *ChannelNumbers

    load sync_times
    load spike_time
    global synct

    [len,widwcv] = size(CSC1_Samples);
    [widlfp] = length(CSC8_Samples(1,:));

    % Get the data.

    wcv = reshape(CSC1_Samples,len*widwcv,1);
    wcv = wcv(synct1id);
    wcv = -wcv; % this sign flip ensures that spikes go upwards.
    wcv = detrend(wcv); % detrend wcv.
    clear CSC1_Samples;

    lf8 = reshape(CSC8_Samples,len*widlfp,1);
    lf8 = lf8(synct2id);
    lf8 = detrend(lf8);

    clear CSC8_Samples;

    synct = synct*10^-6;
    dt = median(diff(synct));

    wcv_minus_spike = wcv;

    %spike subtraction jmm
    spkRmWidth = 20; %width of spike to remove
    for i = 1:length(spkid)
        begPt = wcv(spkid(i)-1);
        endPt = wcv(spkid(i)+spkRmWidth+1);
        interpSlope = (endPt - begPt)/spkRmWidth;
        wcv_minus_spike(spkid(i):spkid(i)+spkRmWidth-1) = [1:spkRmWidth]*interpSlope+begPt;

    end

    %zscore
    wcv_z = wcv_minus_spike - mean(wcv_minus_spike);
    lf8_z = lf8 - mean(lf8);
    wcv_z = wcv_z/std(wcv_z);
    lf8_z = lf8_z/std(lf8_z);

    %filter data
    Fs = median(CSC8_SampleFrequencies);
    lof = .1; hif = 0.5; nyqf = Fs/2;
    lof = lof/nyqf; hif = hif/nyqf;
    [b,a] = butter(2, [lof hif]);
    loflf8 = filtfilt(b,a,lf8_z);
    lofwcv = filtfilt(b,a,wcv_z);

    %higher range bandpass
    lof = .1; hif = 10; nyqf = Fs/2;
    lof = lof/nyqf; hif = hif/nyqf;
    [b,a] = butter(2, [lof hif]);
    medlf8 = filtfilt(b,a,lf8_z);
    medwcv = filtfilt(b,a,wcv_z);


    %take derivative of filtered signals
    diff_loflf8 = [0;diff(loflf8)];
    diff_lofwcv = [0;diff(lofwcv)];
    diff_medlf8 = [0;diff(medlf8)];
    diff_medwcv = [0;diff(medwcv)];

    %zscore derivative signal
    diff_loflf8 = diff_loflf8-mean(diff_loflf8);
    diff_loflf8 = diff_loflf8/std(diff_loflf8);
    diff_lofwcv = diff_lofwcv - mean(diff_lofwcv);
    diff_lofwcv = diff_lofwcv/std(diff_lofwcv);

    diff_medlf8 = diff_medlf8 - mean(diff_medlf8);
    diff_medlf8 = diff_medlf8/std(diff_medlf8);
    diff_medwcv = diff_medwcv - mean(diff_medwcv);
    diff_medwcv = diff_medwcv/std(diff_medwcv);


    %downsample and find peaks
    dsf = 20;
    down_diff_loflf8 = downsample(diff_loflf8,dsf);
    down_diff_lofwcv = downsample(diff_lofwcv,dsf);

    down_diff_medlf8 = downsample(diff_medlf8,dsf);
    down_diff_medwcv = downsample(diff_medwcv,dsf);

    [peakslope_up_wcv,maxslope_up_wcv] = findpeaks(down_diff_lofwcv,'minpeakheight',1,'minpeakdistance',round(0.2/dt/dsf));
    % [peakslope_up_wcv,maxslope_down_wcv] = findpeaks(-down_diff_lofwcv,'minpeakheight',1,'minpeakdistance',round(0.1/dt/dsf));
    [peakslope_up_lf8,maxslope_up_lf8] = findpeaks(down_diff_loflf8,'minpeakheight',1,'minpeakdistance',round(0.2/dt/dsf));


    % go through all up slope peaks
    max_dist = round(2/dt/dsf);

    wcv_up_trans_points = [];

    for i = 1:length(maxslope_up_wcv)

        %find preceding point where slope was below some threshold
        peakthresh = down_diff_medwcv(maxslope_up_wcv(i))/2;
        up_start_temp = find(down_diff_medwcv(1:maxslope_up_wcv(i))<peakthresh,1,'last');

        %%make sure that the signal is below mean at this point
        if medwcv(up_start_temp*dsf) > 0.5
            up_start_temp = [];

            %make sure that this point is within some distance of the slope peak
        elseif abs(up_start_temp - maxslope_up_wcv(i)) > max_dist
            up_start_temp = [];
        end

        wcv_up_trans_points = [wcv_up_trans_points up_start_temp*dsf];

    end


    lf8_up_trans_points = [];

    for i = 1:length(maxslope_up_lf8)

        %find preceding point where slope was below some threshold
        peakthresh = down_diff_medlf8(maxslope_up_lf8(i))/2;
        up_start_temp = find(down_diff_medlf8(1:maxslope_up_lf8(i))<peakthresh,1,'last');

        %%make sure that the signal is below mean at this point
        if medlf8(up_start_temp*dsf) > 0.5
            up_start_temp = [];

            %make sure that this point is within some distance of the slope peak
        elseif abs(up_start_temp - maxslope_up_lf8(i)) > max_dist
            up_start_temp = [];
        end

        lf8_up_trans_points = [lf8_up_trans_points up_start_temp*dsf];

    end


    %get periods
    lf8_period = diff(lf8_up_trans_points)*dt;
    wcv_period = diff(wcv_up_trans_points)*dt;

    %get transition times
    lf8_up_times = synct(lf8_up_trans_points);
    wcv_up_times = synct(wcv_up_trans_points);
    
%     figure
%     set(gcf,'PaperUnits','centimeters');
%     set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
%     set(gcf,'PaperPosition',[0,0,(get(gcf,'PaperSize'))])
%     plot(synct(lf8_up_trans_points(2:end)),fgsmooth(lf8_period,3),'linewidth',2)
%     hold on
%     plot(synct(wcv_up_trans_points(2:end)),fgsmooth(wcv_period,3),'r','linewidth',2)
%     legend('LFP','WCV')
%     xlabel('Time (s)','FontSize',12)
%     ylabel('Local Period (s)','FontSize',12)
%     title('Period vs Time','FontSize',14)
%     tname = ['./jmm_new/period_v_time_' f_names{d}];
%     print('-dpng',tname)
%     tname = ['E:\WC_Germany\JMM_Analysis\period_v_time\' f_names{d}];
%     print('-dpng',tname)
%     close all

    lf8_period_overall = [lf8_period_overall lf8_period];
    wcv_period_overall = [wcv_period_overall wcv_period];
    
end