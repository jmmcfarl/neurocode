
clear all

load E:\WC_Germany\JMM_Analysis\file_directories_jmm.mat
% load E:\WC_Germany\JMM_Analysis\overall_data_jmm.mat

%initialize overall variables
lf8_dci_overall = cell(1,length(dir_array));
lf2_dci_overall = cell(1,length(dir_array));
lf5_dci_overall = cell(1,length(dir_array));
wcv_dci_overall = cell(1,length(dir_array));
%     diff_dci_overall = cell(1,length(dir_array));
corr2_amp_overall = cell(1,length(dir_array));
corr5_amp_overall = cell(1,length(dir_array));
corr8_amp_overall = cell(1,length(dir_array));
corr2_off_overall = cell(1,length(dir_array));
corr5_off_overall = cell(1,length(dir_array));
corr8_off_overall = cell(1,length(dir_array));
corr2_coeff_overall = cell(1,length(dir_array));
corr5_coeff_overall = cell(1,length(dir_array));
corr8_coeff_overall = cell(1,length(dir_array));
lf2_period_overall = cell(1,length(dir_array));
lf5_period_overall = cell(1,length(dir_array));
lf8_period_overall = cell(1,length(dir_array));
wcv_period_overall = cell(1,length(dir_array));
% 
    d = 1;

%cycle through directory tree
while d <= length(dir_array)
    
    cd(dir_array{d})
    pwd

%if this is split data, resplit it
% if dir_array{d}(end-4:end) == 'part1'
%     
%     cd ..
%     save E:\WC_Germany\JMM_Analysis\overall_data_jmm_multLFP.mat *overall d
%     clear d *overall
%     split_very_large_data
%     cd part1
%     load E:\WC_Germany\JMM_Analysis\overall_data_jmm_multLFP.mat *overall d
% 
% end
    
%     if ~exist('./used_data_multLFP.mat')
        
        save E:\WC_Germany\JMM_Analysis\overall_data_jmm_multLFP.mat *overall d
        clear d *overall
        %calculate data
        hip_wc_lfp_spk_shift_jmm
        load E:\WC_Germany\JMM_Analysis\overall_data_jmm_multLFP.mat *overall d
% 
%     else
%         load used_data_multLFP
%     end



    %make sure there is a folder for new figures
    if ~exist('./jmm_8_5')
        !mkdir jmm_8_5
    end


    %zscore
    wcv_z = wcv_minus_spike - mean(wcv_minus_spike);
    lf8_z = lf8 - mean(lf8);
    lf5_z = lf5 - mean(lf5);
    lf2_z = lf2 - mean(lf2);
    wcv_z = wcv_z/std(wcv_z);
    lf8_z = lf8_z/std(lf8_z);
    lf5_z = lf5_z/std(lf5_z);
    lf2_z = lf2_z/std(lf2_z);



    windowSize = 20; %window size in s
    windowSizeBins = round(windowSize/dt);
    timeSlide = 2; %amount to slide time window in s
    binSlide = round(timeSlide/dt);
    numWins = round((length(lf8)-windowSizeBins)/binSlide); %number of windows
    maxLag = round(1/dt); %max lag for cross cov calc
    Fs = mean(CSC8_SampleFrequencies);
    %initializations
    lf8_dci = zeros(1,numWins);
    lf5_dci = zeros(1,numWins);
    lf2_dci = zeros(1,numWins);
    wcv_dci = lf8_dci;
    corr2Amp = lf8_dci;
    corr5Amp = lf8_dci;
    corr8Amp = lf8_dci;
    corr2Off = lf8_dci;
    corr5Off = lf8_dci;
    corr8Off = lf8_dci;
    corr2CoSim = lf8_dci;
    corr5CoSim = lf8_dci;
    corr8CoSim = lf8_dci;

    for i = 1:numWins

        windowCent = (i-1)*binSlide+round(windowSizeBins/2);
        windowBeg = windowCent - round(windowSizeBins/2)+1;
        windowEnd = windowCent + round(windowSizeBins/2);
        lf8Seg = lf8_z(windowBeg:windowEnd);
        lf5Seg = lf5_z(windowBeg:windowEnd);
        lf2Seg = lf2_z(windowBeg:windowEnd);
        %     wcvSeg = wcv(windowBeg:windowEnd);
        wcvSeg_MS = wcv_z(windowBeg:windowEnd);

        lof = .1; hif = 2; nyqf = Fs/2; %set bandpass range
        lof = lof/nyqf; hif = hif/nyqf;
        [b,a] = butter(2, [lof hif]);

        %filter data
        lf8Seg = filtfilt(b,a,lf8Seg);
        lf5Seg = filtfilt(b,a,lf5Seg);
        lf2Seg = filtfilt(b,a,lf2Seg);
        wcvSeg_MS = filtfilt(b,a,wcvSeg_MS);

        %calculate duty cycle index
        lf8_up_time = length(find(lf8Seg>0));
        lf8_down_time = length(find(lf8Seg<0));
        lf5_up_time = length(find(lf5Seg>0));
        lf5_down_time = length(find(lf5Seg<0));
        lf2_up_time = length(find(lf2Seg>0));
        lf2_down_time = length(find(lf2Seg<0));
        wcv_up_time = length(find(wcvSeg_MS>0));
        wcv_down_time = length(find(wcvSeg_MS<-0));

        lf8_dci(i) = (lf8_up_time-lf8_down_time)/(lf8_up_time+lf8_down_time);
        lf5_dci(i) = (lf5_up_time-lf5_down_time)/(lf5_up_time+lf5_down_time);
        lf2_dci(i) = (lf2_up_time-lf2_down_time)/(lf2_up_time+lf2_down_time);

        wcv_dci(i) = (wcv_up_time-wcv_down_time)/(wcv_up_time+wcv_down_time);

%         dci_ratio(i) = wcv_dci(i) - lfp_dci(i);

        %calculate covariance peak and offset

        [corrVec,corrT] = xcov(lf8Seg,wcvSeg_MS,maxLag,'coef');
        [peakAmp,peakLoc] = max(corrVec);
        [peakAbsAmp,peakAbsLoc] = max(abs(corrVec));
        corr8Amp(i) = peakAmp;
        corr8Off(i) = (peakLoc - maxLag)*dt;
%         corrAbsAmp(i) = peakAbsAmp;
%         corrAbsOff(i) = (peakAbsLoc - maxLag)*dt;

        temp = corrcoef(lf8Seg,wcvSeg_MS);
        if ~isempty(temp)
            corr8CoSim(i) = temp(2,1);
        else
            corr8CoSim(i) = NaN;
        end

        [corrVec,corrT] = xcov(lf5Seg,wcvSeg_MS,maxLag,'coef');
        [peakAmp,peakLoc] = max(corrVec);
        [peakAbsAmp,peakAbsLoc] = max(abs(corrVec));
        corr5Amp(i) = peakAmp;
        corr5Off(i) = (peakLoc - maxLag)*dt;
        temp = corrcoef(lf5Seg,wcvSeg_MS);
        if ~isempty(temp)
            corr5CoSim(i) = temp(2,1);
        else
            corr5CoSim(i) = NaN;
        end

        [corrVec,corrT] = xcov(lf2Seg,wcvSeg_MS,maxLag,'coef');
        [peakAmp,peakLoc] = max(corrVec);
        [peakAbsAmp,peakAbsLoc] = max(abs(corrVec));
        corr2Amp(i) = peakAmp;
        corr2Off(i) = (peakLoc - maxLag)*dt;
        temp = corrcoef(lf2Seg,wcvSeg_MS);
        if ~isempty(temp)
            corr2CoSim(i) = temp(2,1);
        else
            corr2CoSim(i) = NaN;
        end

        %

    end



    timeAxis = [1:numWins]*timeSlide;


    %NOW calculate period data

    %filter data
    Fs = median(CSC8_SampleFrequencies);
    lof = .1; hif = 0.5; nyqf = Fs/2;
    lof = lof/nyqf; hif = hif/nyqf;
    [b,a] = butter(2, [lof hif]);
    loflf8 = filtfilt(b,a,lf8_z);
    loflf5 = filtfilt(b,a,lf5_z);
    loflf2 = filtfilt(b,a,lf2_z);
    lofwcv = filtfilt(b,a,wcv_z);

    %higher range bandpass
    lof = .1; hif = 10; nyqf = Fs/2;
    lof = lof/nyqf; hif = hif/nyqf;
    [b,a] = butter(2, [lof hif]);
    medlf8 = filtfilt(b,a,lf8_z);
    medlf5 = filtfilt(b,a,lf5_z);
    medlf2 = filtfilt(b,a,lf2_z);
    medwcv = filtfilt(b,a,wcv_z);


    %take derivative of filtered signals
    diff_loflf8 = [0;diff(loflf8)];
    diff_loflf5 = [0;diff(loflf5)];
    diff_loflf2 = [0;diff(loflf2)];
    diff_lofwcv = [0;diff(lofwcv)];
    diff_medlf8 = [0;diff(medlf8)];
    diff_medlf5 = [0;diff(medlf5)];
    diff_medlf2 = [0;diff(medlf2)];
    diff_medwcv = [0;diff(medwcv)];

    %zscore derivative signal
    diff_loflf8 = diff_loflf8-mean(diff_loflf8);
    diff_loflf8 = diff_loflf8/std(diff_loflf8);
    diff_loflf5 = diff_loflf5-mean(diff_loflf5);
    diff_loflf5 = diff_loflf5/std(diff_loflf5);
    diff_loflf2 = diff_loflf2-mean(diff_loflf2);
    diff_loflf2 = diff_loflf2/std(diff_loflf2);
    diff_lofwcv = diff_lofwcv - mean(diff_lofwcv);
    diff_lofwcv = diff_lofwcv/std(diff_lofwcv);

    diff_medlf8 = diff_medlf8 - mean(diff_medlf8);
    diff_medlf8 = diff_medlf8/std(diff_medlf8);
    diff_medlf5 = diff_medlf5 - mean(diff_medlf5);
    diff_medlf5 = diff_medlf5/std(diff_medlf5);
    diff_medlf2 = diff_medlf2 - mean(diff_medlf2);
    diff_medlf2 = diff_medlf2/std(diff_medlf2);
    diff_medwcv = diff_medwcv - mean(diff_medwcv);
    diff_medwcv = diff_medwcv/std(diff_medwcv);


    %downsample and find peaks
    dsf = 20;
    down_diff_loflf8 = downsample(diff_loflf8,dsf);
    down_diff_loflf5 = downsample(diff_loflf5,dsf);
    down_diff_loflf2 = downsample(diff_loflf2,dsf);
    down_diff_lofwcv = downsample(diff_lofwcv,dsf);

    down_diff_medlf8 = downsample(diff_medlf8,dsf);
    down_diff_medlf5 = downsample(diff_medlf5,dsf);
    down_diff_medlf2 = downsample(diff_medlf2,dsf);
    down_diff_medwcv = downsample(diff_medwcv,dsf);

    [peakslope_up_wcv,maxslope_up_wcv] = findpeaks(down_diff_lofwcv,'minpeakheight',1,'minpeakdistance',round(0.2/dt/dsf));
    % [peakslope_up_wcv,maxslope_down_wcv] = findpeaks(-down_diff_lofwcv,'minpeakheight',1,'minpeakdistance',round(0.1/dt/dsf));
    [peakslope_up_lf8,maxslope_up_lf8] = findpeaks(down_diff_loflf8,'minpeakheight',1,'minpeakdistance',round(0.2/dt/dsf));
    [peakslope_up_lf5,maxslope_up_lf5] = findpeaks(down_diff_loflf5,'minpeakheight',1,'minpeakdistance',round(0.2/dt/dsf));
    [peakslope_up_lf2,maxslope_up_lf2] = findpeaks(down_diff_loflf2,'minpeakheight',1,'minpeakdistance',round(0.2/dt/dsf));

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

    lf5_up_trans_points = [];

    for i = 1:length(maxslope_up_lf5)

        %find preceding point where slope was below some threshold
        peakthresh = down_diff_medlf5(maxslope_up_lf5(i))/2;
        up_start_temp = find(down_diff_medlf5(1:maxslope_up_lf5(i))<peakthresh,1,'last');

        %%make sure that the signal is below mean at this point
        if medlf5(up_start_temp*dsf) > 0.5
            up_start_temp = [];

            %make sure that this point is within some distance of the slope peak
        elseif abs(up_start_temp - maxslope_up_lf5(i)) > max_dist
            up_start_temp = [];
        end

        lf5_up_trans_points = [lf5_up_trans_points up_start_temp*dsf];

    end

    
        lf2_up_trans_points = [];

    for i = 1:length(maxslope_up_lf2)

        %find preceding point where slope was below some threshold
        peakthresh = down_diff_medlf2(maxslope_up_lf2(i))/2;
        up_start_temp = find(down_diff_medlf2(1:maxslope_up_lf2(i))<peakthresh,1,'last');

        %%make sure that the signal is below mean at this point
        if medlf2(up_start_temp*dsf) > 0.5
            up_start_temp = [];

            %make sure that this point is within some distance of the slope peak
        elseif abs(up_start_temp - maxslope_up_lf2(i)) > max_dist
            up_start_temp = [];
        end

        lf2_up_trans_points = [lf2_up_trans_points up_start_temp*dsf];

    end

    
    %get periods
    lf8_period = diff(lf8_up_trans_points)*dt;
    lf5_period = diff(lf5_up_trans_points)*dt;
    lf2_period = diff(lf2_up_trans_points)*dt;
    wcv_period = diff(wcv_up_trans_points)*dt;

    %get transition times
    lf8_up_times = synct(lf8_up_trans_points);
    lf5_up_times = synct(lf5_up_trans_points);
    lf2_up_times = synct(lf2_up_trans_points);
    wcv_up_times = synct(wcv_up_trans_points);


    lf8_dci_overall{d} = lf8_dci;
    lf5_dci_overall{d} = lf5_dci;
    lf2_dci_overall{d} = lf2_dci;
    wcv_dci_overall{d} = wcv_dci;
%     diff_dci_overall{d} = dci_ratio;
    corr8_amp_overall{d} = corr8Amp;
    corr5_amp_overall{d} = corr5Amp;
    corr2_amp_overall{d} = corr2Amp;
    corr8_off_overall{d} = corr8Off;
    corr5_off_overall{d} = corr5Off;
    corr2_off_overall{d} = corr2Off;
    corr8_coeff_overall{d} = corr8CoSim;
    corr5_coeff_overall{d} = corr5CoSim;
    corr2_coeff_overall{d} = corr2CoSim;
    lf8_period_overall{d} = lf8_period;
    lf5_period_overall{d} = lf5_period;
    lf2_period_overall{d} = lf2_period;
    wcv_period_overall{d} = wcv_period;

   
    
% 
%     figure
%     set(gcf,'PaperUnits','centimeters');
%     set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
%     set(gcf,'PaperPosition',[0,0,(get(gcf,'PaperSize'))])
%     plot(synct(lf8_up_trans_points(2:end)),fgsmooth(lf8_period,3),'linewidth',2)
%     hold on
%     plot(synct(lf5_up_trans_points(2:end)),fgsmooth(lf5_period,3),'r','linewidth',2)
%     plot(synct(lf2_up_trans_points(2:end)),fgsmooth(lf2_period,3),'g','linewidth',2)
%     plot(synct(wcv_up_trans_points(2:end)),fgsmooth(wcv_period,3),'k','linewidth',2)
%     legend('LF8','LF5','LF2','WCV')
%     xlabel('Time (s)','FontSize',12)
%     ylabel('Local Period (s)','FontSize',12)
%     title('Period vs Time','FontSize',14)
%     tname = ['./jmm_8_5/period_v_time_' f_names{d}];
%     print('-dpng',tname)
%     tname = ['E:\WC_Germany\JMM_Analysis\multLFP\period_v_time\' f_names{d}];
%     print('-dpng',tname)
%     close all


%     figure
%     set(gcf,'PaperUnits','centimeters');
%     set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
%     set(gcf,'PaperPosition',[0,0,(get(gcf,'PaperSize'))])
%     plot(timeAxis,lf8_dci,'linewidth',2)
%     hold on
%     plot(timeAxis,lf5_dci,'r','linewidth',2)
%     plot(timeAxis,lf2_dci,'g','linewidth',2)
%     plot(timeAxis,wcv_dci,'k','linewidth',2)
%     legend('LF8 DCI','LF5 DCI','LF2 DCI','WCV DCI')
%     xlabel('Time (s)','FontSize',14);
%     ylabel('Duty Cycle Index','FontSize',14)
%     title('Duty Cycle Index Vs Time ','FontSize',16)
%     tname = ['./jmm_8_5/DCI_ratio_sub_' f_names{d}];
%     print('-dpng',tname)
%     tname = ['E:\WC_Germany\JMM_Analysis\multLFP\duty_cycle_index_sub\' f_names{d}];
%     print('-dpng',tname)
%     close all


%     figure
%     set(gcf,'PaperUnits','centimeters');
%     set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
%     set(gcf,'PaperPosition',[0,0,(get(gcf,'PaperSize'))])
%     subplot(3,1,1)
%     plot(timeAxis,corr8Amp,'linewidth',2)
%     hold on
%     plot(timeAxis,corr5Amp,'r','linewidth',2)
%     plot(timeAxis,corr2Amp,'g','linewidth',2)
%     legend('Corr8','Corr5','Corr2')
%     xlabel('Time (s)','FontSize',12)
%     ylabel('Max XCorr','FontSize',12)
%     title('Correlation Peak Amp vs Time','FontSize',14)
%     subplot(3,1,2)    
%     plot(timeAxis,corr8CoSim,'linewidth',2)
%     hold on
%     plot(timeAxis,corr5CoSim,'r','linewidth',2)
%     plot(timeAxis,corr2CoSim,'g','linewidth',2)
%     legend('Corr8','Corr5','Corr2')
%     xlabel('Time (s)','FontSize',12)
%     ylabel('Correlation Coefficient','FontSize',12)
%     title('Correlation Coefficient vs Time','FontSize',14)
%     subplot(3,1,3)
%     plot(timeAxis,corr8Off,'linewidth',2)
%     hold on
%     plot(timeAxis,corr5Off,'r','linewidth',2)
%     plot(timeAxis,corr2Off,'g','linewidth',2)
%     xlabel('Time (s)','FontSize',12)
%     ylabel('Peak Correlation Offset','FontSize',12)
%     title('Correlation Peak Offset vs Time','FontSize',14)
%     tname = ['./jmm_8_5/Correlation_xandco_' f_names{d}];
%     print('-dpng',tname)
%     tname = ['E:\WC_Germany\JMM_Analysis\multLFP\corr_xandco\' f_names{d}];
%     print('-dpng',tname)
%     close all

    
    
    d = d +1

    save E:\WC_Germany\JMM_Analysis\overall_data_jmm_multLFP.mat *overall d

    %prepare for next iteration
    clear all

    load E:\WC_Germany\JMM_Analysis\file_directories_jmm.mat
    load E:\WC_Germany\JMM_Analysis\overall_data_jmm_multLFP.mat

%     d
    
end


cd 'E:\WC_Germany\JMM_Analysis'

%calculate median values of distributions
for i = 1:43
    med_period(i) = median(wcv_period_overall{i});
    med_corr8_amp(i) = median(corr8_amp_overall{i});
    med_corr5_amp(i) = median(corr5_amp_overall{i});
    med_corr2_amp(i) = median(corr2_amp_overall{i});
    med_corr8_coef(i) = median(corr8_coeff_overall{i});
    med_corr5_coef(i) = median(corr5_coeff_overall{i});
    med_corr2_coef(i) = median(corr2_coeff_overall{i});
    med_corr8_off(i) = median(corr8_off_overall{i});
    med_corr5_off(i) = median(corr5_off_overall{i});
    med_corr2_off(i) = median(corr2_off_overall{i});
    med_wcv_dci(i) = median(wcv_dci_overall{i});
end

%initialize total distribution vectors
corr8_amp_total_stell = [];
corr5_amp_total_stell = [];
corr2_amp_total_stell = [];
corr8_coeff_total_stell = [];
corr5_coeff_total_stell = [];
corr2_coeff_total_stell = [];
corr8_off_total_stell = [];
corr5_off_total_stell = [];
corr2_off_total_stell = [];
lf8_dci_total_stell = [];
lf5_dci_total_stell = [];
lf2_dci_total_stell = [];
wcv_dci_total_stell = [];
lf8_period_total_stell = [];
lf5_period_total_stell = [];
lf2_period_total_stell = [];
wcv_period_total_stell = [];
corr8_amp_total_pyr = [];
corr5_amp_total_pyr = [];
corr2_amp_total_pyr = [];
corr8_coeff_total_pyr = [];
corr5_coeff_total_pyr = [];
corr2_coeff_total_pyr = [];
corr8_off_total_pyr = [];
corr5_off_total_pyr = [];
corr2_off_total_pyr = [];
lf8_dci_total_pyr = [];
lf5_dci_total_pyr = [];
lf2_dci_total_pyr = [];
wcv_period_total_pyr = [];
lf8_period_total_pyr = [];
lf5_period_total_pyr = [];
lf2_period_total_pyr = [];
wcv_dci_total_pyr = [];

%for stellate
for i = 1:17
    
    corr8_amp_total_stell = [corr8_amp_total_stell corr8_amp_overall{i}];
    corr5_amp_total_stell = [corr5_amp_total_stell corr5_amp_overall{i}];
    corr2_amp_total_stell = [corr2_amp_total_stell corr2_amp_overall{i}];
    corr8_coeff_total_stell = [corr8_coeff_total_stell corr8_coeff_overall{i}];
    corr5_coeff_total_stell = [corr5_coeff_total_stell corr5_coeff_overall{i}];
    corr2_coeff_total_stell = [corr2_coeff_total_stell corr2_coeff_overall{i}];
    corr8_off_total_stell = [corr8_off_total_stell corr8_off_overall{i}];
    corr5_off_total_stell = [corr5_off_total_stell corr5_off_overall{i}];
    corr2_off_total_stell = [corr2_off_total_stell corr2_off_overall{i}];
    lf8_dci_total_stell = [lf8_dci_total_stell lf8_dci_overall{i}];
    lf5_dci_total_stell = [lf5_dci_total_stell lf5_dci_overall{i}];
    lf2_dci_total_stell = [lf2_dci_total_stell lf2_dci_overall{i}];
    wcv_dci_total_stell = [wcv_dci_total_stell wcv_dci_overall{i}];
    lf8_period_total_stell = [lf8_period_total_stell lf8_period_overall{i}];
    lf5_period_total_stell = [lf5_period_total_stell lf5_period_overall{i}];
    lf2_period_total_stell = [lf2_period_total_stell lf2_period_overall{i}];
    wcv_period_total_stell = [wcv_period_total_stell wcv_period_overall{i}];
    
end


%for pyr
for i = 18:42
    
    corr8_amp_total_pyr = [corr8_amp_total_pyr corr8_amp_overall{i}];
    corr5_amp_total_pyr = [corr5_amp_total_pyr corr5_amp_overall{i}];
    corr2_amp_total_pyr = [corr2_amp_total_pyr corr2_amp_overall{i}];
    corr8_coeff_total_pyr = [corr8_coeff_total_pyr corr8_coeff_overall{i}];
    corr5_coeff_total_pyr = [corr5_coeff_total_pyr corr5_coeff_overall{i}];
    corr2_coeff_total_pyr = [corr2_coeff_total_pyr corr2_coeff_overall{i}];
    corr8_off_total_pyr = [corr8_off_total_pyr corr8_off_overall{i}];
    corr5_off_total_pyr = [corr5_off_total_pyr corr5_off_overall{i}];
    corr2_off_total_pyr = [corr2_off_total_pyr corr2_off_overall{i}];
    lf8_dci_total_pyr = [lf8_dci_total_pyr lf8_dci_overall{i}];
    lf5_dci_total_pyr = [lf5_dci_total_pyr lf5_dci_overall{i}];
    lf2_dci_total_pyr = [lf2_dci_total_pyr lf2_dci_overall{i}];
    wcv_dci_total_pyr = [wcv_dci_total_pyr wcv_dci_overall{i}];
    lf8_period_total_pyr = [lf8_period_total_pyr lf8_period_overall{i}];
    lf5_period_total_pyr = [lf5_period_total_pyr lf5_period_overall{i}];
    lf2_period_total_pyr = [lf2_period_total_pyr lf2_period_overall{i}];
    wcv_period_total_pyr = [wcv_period_total_pyr wcv_period_overall{i}];
    
end

%hist correlation parameters
corr_grid = [-1:0.01:1];
pyr_corr8_amp_hist = hist(corr8_amp_total_pyr,corr_grid);
pyr_corr5_amp_hist = hist(corr5_amp_total_pyr,corr_grid);
pyr_corr2_amp_hist = hist(corr2_amp_total_pyr,corr_grid);

pyr_corr8_coeff_hist = hist(corr8_coeff_total_pyr,corr_grid);
pyr_corr5_coeff_hist = hist(corr5_coeff_total_pyr,corr_grid);
pyr_corr2_coeff_hist = hist(corr2_coeff_total_pyr,corr_grid);

stell_corr8_amp_hist = hist(corr8_amp_total_stell,corr_grid);
stell_corr5_amp_hist = hist(corr5_amp_total_stell,corr_grid);
stell_corr2_amp_hist = hist(corr2_amp_total_stell,corr_grid);

stell_corr8_coeff_hist = hist(corr8_coeff_total_stell,corr_grid);
stell_corr5_coeff_hist = hist(corr5_coeff_total_stell,corr_grid);
stell_corr2_coeff_hist = hist(corr2_coeff_total_stell,corr_grid);

pyr_corr8_amp_hist = pyr_corr8_amp_hist/sum(pyr_corr8_amp_hist);
pyr_corr5_amp_hist = pyr_corr5_amp_hist/sum(pyr_corr5_amp_hist);
pyr_corr2_amp_hist = pyr_corr2_amp_hist/sum(pyr_corr2_amp_hist);

pyr_corr8_coeff_hist = pyr_corr8_coeff_hist/sum(pyr_corr8_coeff_hist);
pyr_corr5_coeff_hist = pyr_corr5_coeff_hist/sum(pyr_corr5_coeff_hist);
pyr_corr2_coeff_hist = pyr_corr2_coeff_hist/sum(pyr_corr2_coeff_hist);

stell_corr8_amp_hist = stell_corr8_amp_hist/sum(stell_corr8_amp_hist);
stell_corr5_amp_hist = stell_corr5_amp_hist/sum(stell_corr5_amp_hist);
stell_corr2_amp_hist = stell_corr2_amp_hist/sum(stell_corr2_amp_hist);

stell_corr8_coeff_hist = stell_corr8_coeff_hist/sum(stell_corr8_coeff_hist);
stell_corr5_coeff_hist = stell_corr5_coeff_hist/sum(stell_corr5_coeff_hist);
stell_corr2_coeff_hist = stell_corr2_coeff_hist/sum(stell_corr2_coeff_hist);

pyr_corr8_off_hist = hist(corr8_off_total_pyr,corr_grid);
pyr_corr5_off_hist = hist(corr5_off_total_pyr,corr_grid);
pyr_corr2_off_hist = hist(corr2_off_total_pyr,corr_grid);

pyr_corr8_off_hist = pyr_corr8_off_hist/sum(pyr_corr8_off_hist);
pyr_corr5_off_hist = pyr_corr5_off_hist/sum(pyr_corr5_off_hist);
pyr_corr2_off_hist = pyr_corr2_off_hist/sum(pyr_corr2_off_hist);

stell_corr8_off_hist = hist(corr8_off_total_stell,corr_grid);
stell_corr5_off_hist = hist(corr5_off_total_stell,corr_grid);
stell_corr2_off_hist = hist(corr2_off_total_stell,corr_grid);

stell_corr8_off_hist = stell_corr8_off_hist/sum(stell_corr8_off_hist);
stell_corr5_off_hist = stell_corr5_off_hist/sum(stell_corr5_off_hist);
stell_corr2_off_hist = stell_corr2_off_hist/sum(stell_corr2_off_hist);



%hist dci parameters
dci_grid = [-1:.01:1];

pyr_lf8_dci_hist = hist(lf8_dci_total_pyr,dci_grid);
pyr_lf5_dci_hist = hist(lf5_dci_total_pyr,dci_grid);
pyr_lf2_dci_hist = hist(lf2_dci_total_pyr,dci_grid);

pyr_lf8_dci_hist = pyr_lf8_dci_hist/sum(pyr_lf8_dci_hist);
pyr_lf5_dci_hist = pyr_lf5_dci_hist/sum(pyr_lf5_dci_hist);
pyr_lf2_dci_hist = pyr_lf2_dci_hist/sum(pyr_lf2_dci_hist);

pyr_wcv_dci_hist = hist(wcv_dci_total_pyr,dci_grid);
pyr_wcv_dci_hist = pyr_wcv_dci_hist/sum(pyr_wcv_dci_hist);

stell_lf8_dci_hist = hist(lf8_dci_total_stell,dci_grid);
stell_lf5_dci_hist = hist(lf5_dci_total_stell,dci_grid);
stell_lf2_dci_hist = hist(lf2_dci_total_stell,dci_grid);

stell_lf8_dci_hist = stell_lf8_dci_hist/sum(stell_lf8_dci_hist);
stell_lf5_dci_hist = stell_lf5_dci_hist/sum(stell_lf5_dci_hist);
stell_lf2_dci_hist = stell_lf2_dci_hist/sum(stell_lf2_dci_hist);

stell_wcv_dci_hist = hist(wcv_dci_total_stell,dci_grid);
stell_wcv_dci_hist = stell_wcv_dci_hist/sum(stell_wcv_dci_hist);

total_lf8 = [lf8_dci_total_pyr lf8_dci_total_stell];
total_lf5 = [lf5_dci_total_pyr lf5_dci_total_stell];
total_lf2 = [lf2_dci_total_pyr lf2_dci_total_stell];

total_lf8_hist = hist(total_lf8,dci_grid);
total_lf5_hist = hist(total_lf5,dci_grid);
total_lf2_hist = hist(total_lf2,dci_grid);

total_lf8_hist = total_lf8_hist/sum(total_lf8_hist);
total_lf5_hist = total_lf5_hist/sum(total_lf5_hist);
total_lf2_hist = total_lf2_hist/sum(total_lf2_hist);


%hist period params
period_grid = [0:.01:60];
pyr_lf8_period_hist = hist(lf8_period_total_pyr,period_grid);
pyr_lf8_period_hist = pyr_lf8_period_hist/sum(pyr_lf8_period_hist);
pyr_wcv_period_hist = hist(wcv_period_total_pyr,period_grid);
pyr_wcv_period_hist = pyr_wcv_period_hist/sum(pyr_wcv_period_hist);
stell_lf8_period_hist = hist(lf8_period_total_stell,period_grid);
stell_lf8_period_hist = stell_lf8_period_hist/sum(stell_lf8_period_hist);
stell_wcv_period_hist = hist(wcv_period_total_stell,period_grid);
stell_wcv_period_hist = stell_wcv_period_hist/sum(stell_wcv_period_hist);



save overall_data_jmm_multLFP


% figure
% subplot(3,1,1) 
% plot(corr_grid,pyr_corr8_amp_hist,'linewidth',2)
% hold on
% plot(corr_grid,stell_corr8_amp_hist,'r','linewidth',2)
% xlabel('Peak Cross Correlation','FontSize',12)
% ylabel('Probability','FontSize',12)
% legend('Pyramidal','Stellate')
% title('Max XCorr LF8','FontSize',14)
% subplot(3,1,2)
% plot(corr_grid,pyr_corr5_amp_hist,'linewidth',2)
% hold on
% plot(corr_grid,stell_corr5_amp_hist,'r','linewidth',2)
% xlabel('Peak Cross Correlation','FontSize',12)
% ylabel('Probability','FontSize',12)
% legend('Pyramidal','Stellate')
% title('Max XCorr LF5','FontSize',14)
% subplot(3,1,3)
% plot(corr_grid,pyr_corr2_amp_hist,'linewidth',2)
% hold on
% plot(corr_grid,stell_corr2_amp_hist,'r','linewidth',2)
% xlabel('Peak Cross Correlation','FontSize',12)
% ylabel('Probability','FontSize',12)
% legend('Pyramidal','Stellate')
% title('Max XCorr LF2','FontSize',14)
% 
% 
% 
% figure
% subplot(3,1,1)
% plot(corr_grid,pyr_corr8_coeff_hist,'linewidth',2)
% hold on
% plot(corr_grid,stell_corr8_coeff_hist,'r','linewidth',2)
% xlabel('Correlation','FontSize',12)
% ylabel('Probability','FontSize',12)
% legend('Pyramidal','Stellate')
% title('Correlation LF8','FontSize',14)
% subplot(3,1,2)
% plot(corr_grid,pyr_corr5_coeff_hist,'linewidth',2)
% hold on
% plot(corr_grid,stell_corr5_coeff_hist,'r','linewidth',2)
% xlabel('Correlation','FontSize',12)
% ylabel('Probability','FontSize',12)
% legend('Pyramidal','Stellate')
% title('Correlation LF5','FontSize',14)
% subplot(3,1,3)
% plot(corr_grid,pyr_corr2_coeff_hist,'linewidth',2)
% hold on
% plot(corr_grid,stell_corr2_coeff_hist,'r','linewidth',2)
% xlabel('Correlation','FontSize',12)
% ylabel('Probability','FontSize',12)
% legend('Pyramidal','Stellate')
% title('Correlation LF2','FontSize',14)
% 
% figure
% subplot(3,1,1)
% plot(corr_grid,pyr_corr8_off_hist,'linewidth',2)
% hold on
% plot(corr_grid,stell_corr8_off_hist,'r','linewidth',2)
% xlabel('Peak XCorr Offset (s)','FontSize',12)
% ylabel('Probability','FontSize',12)
% legend('Pyramidal','Stellate')
% title('Peak Xcorr Offset LF8','FontSize',14)
% subplot(3,1,2)
% plot(corr_grid,pyr_corr5_off_hist,'linewidth',2)
% hold on
% plot(corr_grid,stell_corr5_off_hist,'r','linewidth',2)
% xlabel('Peak XCorr Offset (s)','FontSize',12)
% ylabel('Probability','FontSize',12)
% legend('Pyramidal','Stellate')
% title('Peak Xcorr Offset LF5','FontSize',14)
% subplot(3,1,3)
% plot(corr_grid,pyr_corr2_off_hist,'linewidth',2)
% hold on
% plot(corr_grid,stell_corr2_off_hist,'r','linewidth',2)
% xlabel('Peak XCorr Offset (s)','FontSize',12)
% ylabel('Probability','FontSize',12)
% legend('Pyramidal','Stellate')
% title('Peak Xcorr Offset LF2','FontSize',14)
% 
% figure
% subplot(3,1,1)
% plot(dci_grid,pyr_lf8_dci_hist,'linewidth',2)
% hold on
% plot(dci_grid,stell_lf8_dci_hist,'r','linewidth',2)
% legend('Pyramidal','Stellate')
% xlabel('DCI')
% ylabel('Probability','FontSize',14)
% title('LF8 DCI Distributions','FontSize',16)
% subplot(3,1,2)
% plot(dci_grid,pyr_lf5_dci_hist,'linewidth',2)
% hold on
% plot(dci_grid,stell_lf5_dci_hist,'r','linewidth',2)
% legend('Pyramidal','Stellate')
% xlabel('DCI')
% ylabel('Probability','FontSize',14)
% title('LF5 DCI Distributions','FontSize',16)
% subplot(3,1,3)
% plot(dci_grid,pyr_lf2_dci_hist,'linewidth',2)
% hold on
% plot(dci_grid,stell_lf2_dci_hist,'r','linewidth',2)
% legend('Pyramidal','Stellate')
% xlabel('DCI')
% ylabel('Probability','FontSize',14)
% title('LF2 DCI Distributions','FontSize',16)
% 
% figure
% plot(dci_grid,pyr_wcv_dci_hist,'linewidth',2)
% hold on
% plot(dci_grid,stell_wcv_dci_hist,'r','linewidth',2)
% legend('Pyramidal','Stellate')
% xlabel('DCI')
% ylabel('Probability','FontSize',14)
% title('WCV DCI Distributions','FontSize',16)
% plot(dci_grid,total_lfp_hist,'k','linewidth',2)
% legend('Pyramidal','Stellate','AVG LFP')
% 
% 
% figure
% plot(period_grid,fgsmooth(pyr_wcv_period_hist,3),'linewidth',2)
% hold on
% plot(period_grid,fgsmooth(stell_wcv_period_hist,3),'r','linewidth',2)
% shg
% plot(period_grid,fgsmooth(stell_lf8_period_hist,3),'k','linewidth',2)
% xlabel('Period (Hz)')
% ylabel('Probability','FontSize',14)
% title('Period Distributions','FontSize',16)
% legend('Pyramidal WCV','Stellate WCV','AVG LFP')
% 
% 
% %check for depth dependance
% load depths
% 
% figure
% plot(depth(1:17),med_period(1:17),'o')
% hold on
% plot(depth(17:end),med_period(17:end),'ro')
% legend('stellate','pyramidal')
% title('Depth v period','FontSize',14)
% 
% figure
% subplot(3,1,1)
% plot(depth(1:17),med_corr8_amp(1:17),'o')
% hold on
% plot(depth(17:end),med_corr8_amp(17:end),'ro')
% legend('stellate','pyramidal')
% title('Depth v peak xcorr LF8','FontSize',14)
% subplot(3,1,2)
% plot(depth(1:17),med_corr5_amp(1:17),'o')
% hold on
% plot(depth(17:end),med_corr5_amp(17:end),'ro')
% legend('stellate','pyramidal')
% title('Depth v peak xcorr LF5','FontSize',14)
% subplot(3,1,3)
% plot(depth(1:17),med_corr2_amp(1:17),'o')
% hold on
% plot(depth(17:end),med_corr2_amp(17:end),'ro')
% legend('stellate','pyramidal')
% title('Depth v peak xcorr LF2','FontSize',14)
% 
% figure
% subplot(3,1,1)
% plot(depth(1:17),med_corr8_coef(1:17),'o')
% hold on
% plot(depth(17:end),med_corr8_coef(17:end),'ro')
% legend('stellate','pyramidal')
% title('Depth v corr coef LF8','FontSize',14)
% subplot(3,1,2)
% plot(depth(1:17),med_corr5_coef(1:17),'o')
% hold on
% plot(depth(17:end),med_corr5_coef(17:end),'ro')
% legend('stellate','pyramidal')
% title('Depth v corr coef LF5','FontSize',14)
% subplot(3,1,3)
% plot(depth(1:17),med_corr2_coef(1:17),'o')
% hold on
% plot(depth(17:end),med_corr2_coef(17:end),'ro')
% legend('stellate','pyramidal')
% title('Depth v corr coef LF2','FontSize',14)
% 
% figure
% subplot(3,1,1)
% plot(depth(1:17),med_corr8_off(1:17),'o')
% hold on
% plot(depth(17:end),med_corr8_off(17:end),'ro')
% legend('stellate','pyramidal')
% title('Depth v corr offset LF8','FontSize',14)
% subplot(3,1,2)
% plot(depth(1:17),med_corr5_off(1:17),'o')
% hold on
% plot(depth(17:end),med_corr5_off(17:end),'ro')
% legend('stellate','pyramidal')
% title('Depth v corr offset LF5','FontSize',14)
% subplot(3,1,3)
% plot(depth(1:17),med_corr2_off(1:17),'o')
% hold on
% plot(depth(17:end),med_corr2_off(17:end),'ro')
% legend('stellate','pyramidal')
% title('Depth v corr offset LF2','FontSize',14)
% 
% figure
% plot(depth(1:17),med_wcv_dci(1:17),'o')
% hold on
% plot(depth(17:end),med_wcv_dci(17:end),'ro')
% legend('stellate','pyramidal')
% title('Depth v DCI','FontSize',14)

