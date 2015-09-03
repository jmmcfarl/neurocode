
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

% % if this is split data, resplit it
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
        
        save E:\WC_Germany\JMM_Analysis\overall_data_jmm_multLFP_2.mat *overall d
        clear d *overall
        %calculate data
        hip_wc_lfp_spk_shift_jmm
        load E:\WC_Germany\JMM_Analysis\overall_data_jmm_multLFP_2.mat *overall d
% 
%     else
%         load used_data_multLFP
%     end



    %make sure there is a folder for new figures
    if ~exist('./jmm_8_6')
        !mkdir jmm_8_6
    end


    %zscore
    wcv_z = wcv_minus_spike - mean(wcv_minus_spike);
    lf3_z = lf3 - mean(lf3);
    lf4_z = lf4 - mean(lf4);
    lf7_z = lf7 - mean(lf7);
    wcv_z = wcv_z/std(wcv_z);
    lf3_z = lf3_z/std(lf3_z);
    lf4_z = lf4_z/std(lf4_z);
    lf7_z = lf7_z/std(lf7_z);



    windowSize = 20; %window size in s
    windowSizeBins = round(windowSize/dt);
    timeSlide = 2; %amount to slide time window in s
    binSlide = round(timeSlide/dt);
    numWins = round((length(lf3)-windowSizeBins)/binSlide); %number of windows
    maxLag = round(1/dt); %max lag for cross cov calc

    if exist('CSC3_SampleFrequencies')
        Fs = mean(CSC3_SampleFrequencies);
    else
        Fs = mean(CSC8_SampleFrequencies);
    end
    
    %initializations
    lf3_dci = zeros(1,numWins);
    lf4_dci = zeros(1,numWins);
    lf7_dci = zeros(1,numWins);
    wcv_dci = lf3_dci;
    corr3Amp = lf3_dci;
    corr4Amp = lf3_dci;
    corr7Amp = lf3_dci;
    corr3Off = lf3_dci;
    corr4Off = lf3_dci;
    corr7Off = lf3_dci;
    corr3CoSim = lf3_dci;
    corr4CoSim = lf3_dci;
    corr7CoSim = lf3_dci;

    for i = 1:numWins

        windowCent = (i-1)*binSlide+round(windowSizeBins/2);
        windowBeg = windowCent - round(windowSizeBins/2)+1;
        windowEnd = windowCent + round(windowSizeBins/2);
        lf3Seg = lf3_z(windowBeg:windowEnd);
        lf4Seg = lf4_z(windowBeg:windowEnd);
        lf7Seg = lf7_z(windowBeg:windowEnd);
        %     wcvSeg = wcv(windowBeg:windowEnd);
        wcvSeg_MS = wcv_z(windowBeg:windowEnd);

        lof = .1; hif = 2; nyqf = Fs/2; %set bandpass range
        lof = lof/nyqf; hif = hif/nyqf;
        [b,a] = butter(2, [lof hif]);

        %filter data
        lf3Seg = filtfilt(b,a,lf3Seg);
        lf4Seg = filtfilt(b,a,lf4Seg);
        lf7Seg = filtfilt(b,a,lf7Seg);
        wcvSeg_MS = filtfilt(b,a,wcvSeg_MS);

        %calculate duty cycle index
        lf3_up_time = length(find(lf3Seg>0));
        lf3_down_time = length(find(lf3Seg<0));
        lf4_up_time = length(find(lf4Seg>0));
        lf4_down_time = length(find(lf4Seg<0));
        lf7_up_time = length(find(lf7Seg>0));
        lf7_down_time = length(find(lf7Seg<0));
        wcv_up_time = length(find(wcvSeg_MS>0));
        wcv_down_time = length(find(wcvSeg_MS<-0));

        lf3_dci(i) = (lf3_up_time-lf3_down_time)/(lf3_up_time+lf3_down_time);
        lf4_dci(i) = (lf4_up_time-lf4_down_time)/(lf4_up_time+lf4_down_time);
        lf7_dci(i) = (lf7_up_time-lf7_down_time)/(lf7_up_time+lf7_down_time);

        wcv_dci(i) = (wcv_up_time-wcv_down_time)/(wcv_up_time+wcv_down_time);

%         dci_ratio(i) = wcv_dci(i) - lfp_dci(i);

        %calculate covariance peak and offset

        [corrVec,corrT] = xcov(lf3Seg,wcvSeg_MS,maxLag,'coef');
        [peakAmp,peakLoc] = max(corrVec);
        [peakAbsAmp,peakAbsLoc] = max(abs(corrVec));
        corr3Amp(i) = peakAmp;
        corr3Off(i) = (peakLoc - maxLag)*dt;
%         corrAbsAmp(i) = peakAbsAmp;
%         corrAbsOff(i) = (peakAbsLoc - maxLag)*dt;

        temp = corrcoef(lf3Seg,wcvSeg_MS);
        if ~isempty(temp)
            corr3CoSim(i) = temp(2,1);
        else
            corr3CoSim(i) = NaN;
        end

        [corrVec,corrT] = xcov(lf4Seg,wcvSeg_MS,maxLag,'coef');
        [peakAmp,peakLoc] = max(corrVec);
        [peakAbsAmp,peakAbsLoc] = max(abs(corrVec));
        corr4Amp(i) = peakAmp;
        corr4Off(i) = (peakLoc - maxLag)*dt;
        temp = corrcoef(lf4Seg,wcvSeg_MS);
        if ~isempty(temp)
            corr4CoSim(i) = temp(2,1);
        else
            corr4CoSim(i) = NaN;
        end

        [corrVec,corrT] = xcov(lf7Seg,wcvSeg_MS,maxLag,'coef');
        [peakAmp,peakLoc] = max(corrVec);
        [peakAbsAmp,peakAbsLoc] = max(abs(corrVec));
        corr7Amp(i) = peakAmp;
        corr7Off(i) = (peakLoc - maxLag)*dt;
        temp = corrcoef(lf7Seg,wcvSeg_MS);
        if ~isempty(temp)
            corr7CoSim(i) = temp(2,1);
        else
            corr7CoSim(i) = NaN;
        end

        %

    end



    timeAxis = [1:numWins]*timeSlide;


    %NOW calculate period data

    %filter data
%     Fs = median(CSC3_SampleFrequencies);
    lof = .1; hif = 0.5; nyqf = Fs/2;
    lof = lof/nyqf; hif = hif/nyqf;
    [b,a] = butter(2, [lof hif]);
    loflf3 = filtfilt(b,a,lf3_z);
    loflf4 = filtfilt(b,a,lf4_z);
    loflf7 = filtfilt(b,a,lf7_z);
    lofwcv = filtfilt(b,a,wcv_z);

    %higher range bandpass
    lof = .1; hif = 10; nyqf = Fs/2;
    lof = lof/nyqf; hif = hif/nyqf;
    [b,a] = butter(2, [lof hif]);
    medlf3 = filtfilt(b,a,lf3_z);
    medlf4 = filtfilt(b,a,lf4_z);
    medlf7 = filtfilt(b,a,lf7_z);
    medwcv = filtfilt(b,a,wcv_z);


    %take derivative of filtered signals
    diff_loflf3 = [0;diff(loflf3)];
    diff_loflf4 = [0;diff(loflf4)];
    diff_loflf7 = [0;diff(loflf7)];
    diff_lofwcv = [0;diff(lofwcv)];
    diff_medlf3 = [0;diff(medlf3)];
    diff_medlf4 = [0;diff(medlf4)];
    diff_medlf7 = [0;diff(medlf7)];
    diff_medwcv = [0;diff(medwcv)];

    %zscore derivative signal
    diff_loflf3 = diff_loflf3-mean(diff_loflf3);
    diff_loflf3 = diff_loflf3/std(diff_loflf3);
    diff_loflf4 = diff_loflf4-mean(diff_loflf4);
    diff_loflf4 = diff_loflf4/std(diff_loflf4);
    diff_loflf7 = diff_loflf7-mean(diff_loflf7);
    diff_loflf7 = diff_loflf7/std(diff_loflf7);
    diff_lofwcv = diff_lofwcv - mean(diff_lofwcv);
    diff_lofwcv = diff_lofwcv/std(diff_lofwcv);

    diff_medlf3 = diff_medlf3 - mean(diff_medlf3);
    diff_medlf3 = diff_medlf3/std(diff_medlf3);
    diff_medlf4 = diff_medlf4 - mean(diff_medlf4);
    diff_medlf4 = diff_medlf4/std(diff_medlf4);
    diff_medlf7 = diff_medlf7 - mean(diff_medlf7);
    diff_medlf7 = diff_medlf7/std(diff_medlf7);
    diff_medwcv = diff_medwcv - mean(diff_medwcv);
    diff_medwcv = diff_medwcv/std(diff_medwcv);


    %downsample and find peaks
    dsf = 20;
    down_diff_loflf3 = downsample(diff_loflf3,dsf);
    down_diff_loflf4 = downsample(diff_loflf4,dsf);
    down_diff_loflf7 = downsample(diff_loflf7,dsf);
    down_diff_lofwcv = downsample(diff_lofwcv,dsf);

    down_diff_medlf3 = downsample(diff_medlf3,dsf);
    down_diff_medlf4 = downsample(diff_medlf4,dsf);
    down_diff_medlf7 = downsample(diff_medlf7,dsf);
    down_diff_medwcv = downsample(diff_medwcv,dsf);

    [peakslope_up_wcv,maxslope_up_wcv] = findpeaks(down_diff_lofwcv,'minpeakheight',1,'minpeakdistance',round(0.2/dt/dsf));
    % [peakslope_up_wcv,maxslope_down_wcv] = findpeaks(-down_diff_lofwcv,'minpeakheight',1,'minpeakdistance',round(0.1/dt/dsf));
    [peakslope_up_lf3,maxslope_up_lf3] = findpeaks(down_diff_loflf3,'minpeakheight',1,'minpeakdistance',round(0.2/dt/dsf));
    [peakslope_up_lf4,maxslope_up_lf4] = findpeaks(down_diff_loflf4,'minpeakheight',1,'minpeakdistance',round(0.2/dt/dsf));
    [peakslope_up_lf7,maxslope_up_lf7] = findpeaks(down_diff_loflf7,'minpeakheight',1,'minpeakdistance',round(0.2/dt/dsf));

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


    lf3_up_trans_points = [];

    for i = 1:length(maxslope_up_lf3)

        %find preceding point where slope was below some threshold
        peakthresh = down_diff_medlf3(maxslope_up_lf3(i))/2;
        up_start_temp = find(down_diff_medlf3(1:maxslope_up_lf3(i))<peakthresh,1,'last');

        %%make sure that the signal is below mean at this point
        if medlf3(up_start_temp*dsf) > 0.5
            up_start_temp = [];

            %make sure that this point is within some distance of the slope peak
        elseif abs(up_start_temp - maxslope_up_lf3(i)) > max_dist
            up_start_temp = [];
        end

        lf3_up_trans_points = [lf3_up_trans_points up_start_temp*dsf];

    end

    lf4_up_trans_points = [];

    for i = 1:length(maxslope_up_lf4)

        %find preceding point where slope was below some threshold
        peakthresh = down_diff_medlf4(maxslope_up_lf4(i))/2;
        up_start_temp = find(down_diff_medlf4(1:maxslope_up_lf4(i))<peakthresh,1,'last');

        %%make sure that the signal is below mean at this point
        if medlf4(up_start_temp*dsf) > 0.5
            up_start_temp = [];

            %make sure that this point is within some distance of the slope peak
        elseif abs(up_start_temp - maxslope_up_lf4(i)) > max_dist
            up_start_temp = [];
        end

        lf4_up_trans_points = [lf4_up_trans_points up_start_temp*dsf];

    end

    
        lf7_up_trans_points = [];

    for i = 1:length(maxslope_up_lf7)

        %find preceding point where slope was below some threshold
        peakthresh = down_diff_medlf7(maxslope_up_lf7(i))/2;
        up_start_temp = find(down_diff_medlf7(1:maxslope_up_lf7(i))<peakthresh,1,'last');

        %%make sure that the signal is below mean at this point
        if medlf7(up_start_temp*dsf) > 0.5
            up_start_temp = [];

            %make sure that this point is within some distance of the slope peak
        elseif abs(up_start_temp - maxslope_up_lf7(i)) > max_dist
            up_start_temp = [];
        end

        lf7_up_trans_points = [lf7_up_trans_points up_start_temp*dsf];

    end

    
    %get periods
    lf3_period = diff(lf3_up_trans_points)*dt;
    lf4_period = diff(lf4_up_trans_points)*dt;
    lf7_period = diff(lf7_up_trans_points)*dt;
    wcv_period = diff(wcv_up_trans_points)*dt;

    %get transition times
    lf3_up_times = synct(lf3_up_trans_points);
    lf4_up_times = synct(lf4_up_trans_points);
    lf7_up_times = synct(lf7_up_trans_points);
    wcv_up_times = synct(wcv_up_trans_points);


    lf3_dci_overall{d} = lf3_dci;
    lf4_dci_overall{d} = lf4_dci;
    lf7_dci_overall{d} = lf7_dci;
    wcv_dci_overall{d} = wcv_dci;
%     diff_dci_overall{d} = dci_ratio;
    corr3_amp_overall{d} = corr3Amp;
    corr4_amp_overall{d} = corr4Amp;
    corr7_amp_overall{d} = corr7Amp;
    corr3_off_overall{d} = corr3Off;
    corr4_off_overall{d} = corr4Off;
    corr7_off_overall{d} = corr7Off;
    corr3_coeff_overall{d} = corr3CoSim;
    corr4_coeff_overall{d} = corr4CoSim;
    corr7_coeff_overall{d} = corr7CoSim;
    lf3_period_overall{d} = lf3_period;
    lf4_period_overall{d} = lf4_period;
    lf7_period_overall{d} = lf7_period;
    wcv_period_overall{d} = wcv_period;

   
    

    figure
    set(gcf,'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
    set(gcf,'PaperPosition',[0,0,(get(gcf,'PaperSize'))])
    plot(synct(lf3_up_trans_points(2:end)),fgsmooth(lf3_period,3),'linewidth',2)
    hold on
    plot(synct(lf4_up_trans_points(2:end)),fgsmooth(lf4_period,3),'r','linewidth',2)
    plot(synct(lf7_up_trans_points(2:end)),fgsmooth(lf7_period,3),'g','linewidth',2)
    plot(synct(wcv_up_trans_points(2:end)),fgsmooth(wcv_period,3),'k','linewidth',2)
    legend('LF3','LF4','LF7','WCV')
    xlabel('Time (s)','FontSize',12)
    ylabel('Local Period (s)','FontSize',12)
    title('Period vs Time','FontSize',14)
    tname = ['./jmm_8_6/period_v_time_' f_names{d}];
    print('-dpng',tname)
    tname = ['E:\WC_Germany\JMM_Analysis\multLFP_2\period_v_time\' f_names{d}];
    print('-dpng',tname)
    close all


    figure
    set(gcf,'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
    set(gcf,'PaperPosition',[0,0,(get(gcf,'PaperSize'))])
    plot(timeAxis,lf3_dci,'linewidth',2)
    hold on
    plot(timeAxis,lf4_dci,'r','linewidth',2)
    plot(timeAxis,lf7_dci,'g','linewidth',2)
    plot(timeAxis,wcv_dci,'k','linewidth',2)
    legend('LF3 DCI','LF4 DCI','LF7 DCI','WCV DCI')
    xlabel('Time (s)','FontSize',14);
    ylabel('Duty Cycle Index','FontSize',14)
    title('Duty Cycle Index Vs Time ','FontSize',16)
    tname = ['./jmm_8_6/DCI_ratio_sub_' f_names{d}];
    print('-dpng',tname)
    tname = ['E:\WC_Germany\JMM_Analysis\multLFP_2\duty_cycle_index_sub\' f_names{d}];
    print('-dpng',tname)
    close all


    figure
    set(gcf,'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
    set(gcf,'PaperPosition',[0,0,(get(gcf,'PaperSize'))])
    subplot(3,1,1)
    plot(timeAxis,corr3Amp,'linewidth',2)
    hold on
    plot(timeAxis,corr4Amp,'r','linewidth',2)
    plot(timeAxis,corr7Amp,'g','linewidth',2)
    legend('Corr3','Corr4','Corr7')
    xlabel('Time (s)','FontSize',12)
    ylabel('Max XCorr','FontSize',12)
    title('Correlation Peak Amp vs Time','FontSize',14)
    subplot(3,1,2)    
    plot(timeAxis,corr3CoSim,'linewidth',2)
    hold on
    plot(timeAxis,corr4CoSim,'r','linewidth',2)
    plot(timeAxis,corr7CoSim,'g','linewidth',2)
    legend('Corr3','Corr4','Corr7')
    xlabel('Time (s)','FontSize',12)
    ylabel('Correlation Coefficient','FontSize',12)
    title('Correlation Coefficient vs Time','FontSize',14)
    subplot(3,1,3)
    plot(timeAxis,corr3Off,'linewidth',2)
    hold on
    plot(timeAxis,corr4Off,'r','linewidth',2)
    plot(timeAxis,corr7Off,'g','linewidth',2)
    xlabel('Time (s)','FontSize',12)
    ylabel('Peak Correlation Offset','FontSize',12)
    title('Correlation Peak Offset vs Time','FontSize',14)
    tname = ['./jmm_8_6/Correlation_xandco_' f_names{d}];
    print('-dpng',tname)
    tname = ['E:\WC_Germany\JMM_Analysis\multLFP_2\corr_xandco\' f_names{d}];
    print('-dpng',tname)
    close all

    
    
    d = d +1

    save E:\WC_Germany\JMM_Analysis\overall_data_jmm_multLFP_2.mat *overall d

    %prepare for next iteration
    clear all

    load E:\WC_Germany\JMM_Analysis\file_directories_jmm.mat
    load E:\WC_Germany\JMM_Analysis\overall_data_jmm_multLFP_2.mat

%     d
    
end


cd 'E:\WC_Germany\JMM_Analysis'

%calculate median values of distributions
for i = 1:43
    med_period(i) = median(wcv_period_overall{i});
    med_corr3_amp(i) = median(corr3_amp_overall{i});
    med_corr4_amp(i) = median(corr4_amp_overall{i});
    med_corr7_amp(i) = median(corr7_amp_overall{i});
    med_corr3_coef(i) = median(corr3_coeff_overall{i});
    med_corr4_coef(i) = median(corr4_coeff_overall{i});
    med_corr7_coef(i) = median(corr7_coeff_overall{i});
    med_corr3_off(i) = median(corr3_off_overall{i});
    med_corr4_off(i) = median(corr4_off_overall{i});
    med_corr7_off(i) = median(corr7_off_overall{i});
    med_wcv_dci(i) = median(wcv_dci_overall{i});
end

%initialize total distribution vectors
corr3_amp_total_stell = [];
corr4_amp_total_stell = [];
corr7_amp_total_stell = [];
corr3_coeff_total_stell = [];
corr4_coeff_total_stell = [];
corr7_coeff_total_stell = [];
corr3_off_total_stell = [];
corr4_off_total_stell = [];
corr7_off_total_stell = [];
lf3_dci_total_stell = [];
lf4_dci_total_stell = [];
lf7_dci_total_stell = [];
wcv_dci_total_stell = [];
lf3_period_total_stell = [];
lf4_period_total_stell = [];
lf7_period_total_stell = [];
wcv_period_total_stell = [];
corr3_amp_total_pyr = [];
corr4_amp_total_pyr = [];
corr7_amp_total_pyr = [];
corr3_coeff_total_pyr = [];
corr4_coeff_total_pyr = [];
corr7_coeff_total_pyr = [];
corr3_off_total_pyr = [];
corr4_off_total_pyr = [];
corr7_off_total_pyr = [];
lf3_dci_total_pyr = [];
lf4_dci_total_pyr = [];
lf7_dci_total_pyr = [];
wcv_period_total_pyr = [];
lf3_period_total_pyr = [];
lf4_period_total_pyr = [];
lf7_period_total_pyr = [];
wcv_dci_total_pyr = [];

%for stellate
for i = 1:17
    
    corr3_amp_total_stell = [corr3_amp_total_stell corr3_amp_overall{i}];
    corr4_amp_total_stell = [corr4_amp_total_stell corr4_amp_overall{i}];
    corr7_amp_total_stell = [corr7_amp_total_stell corr7_amp_overall{i}];
    corr3_coeff_total_stell = [corr3_coeff_total_stell corr3_coeff_overall{i}];
    corr4_coeff_total_stell = [corr4_coeff_total_stell corr4_coeff_overall{i}];
    corr7_coeff_total_stell = [corr7_coeff_total_stell corr7_coeff_overall{i}];
    corr3_off_total_stell = [corr3_off_total_stell corr3_off_overall{i}];
    corr4_off_total_stell = [corr4_off_total_stell corr4_off_overall{i}];
    corr7_off_total_stell = [corr7_off_total_stell corr7_off_overall{i}];
    lf3_dci_total_stell = [lf3_dci_total_stell lf3_dci_overall{i}];
    lf4_dci_total_stell = [lf4_dci_total_stell lf4_dci_overall{i}];
    lf7_dci_total_stell = [lf7_dci_total_stell lf7_dci_overall{i}];
    wcv_dci_total_stell = [wcv_dci_total_stell wcv_dci_overall{i}];
    lf3_period_total_stell = [lf3_period_total_stell lf3_period_overall{i}];
    lf4_period_total_stell = [lf4_period_total_stell lf4_period_overall{i}];
    lf7_period_total_stell = [lf7_period_total_stell lf7_period_overall{i}];
    wcv_period_total_stell = [wcv_period_total_stell wcv_period_overall{i}];
    
end


%for pyr
for i = 18:43
    
    corr3_amp_total_pyr = [corr3_amp_total_pyr corr3_amp_overall{i}];
    corr4_amp_total_pyr = [corr4_amp_total_pyr corr4_amp_overall{i}];
    corr7_amp_total_pyr = [corr7_amp_total_pyr corr7_amp_overall{i}];
    corr3_coeff_total_pyr = [corr3_coeff_total_pyr corr3_coeff_overall{i}];
    corr4_coeff_total_pyr = [corr4_coeff_total_pyr corr4_coeff_overall{i}];
    corr7_coeff_total_pyr = [corr7_coeff_total_pyr corr7_coeff_overall{i}];
    corr3_off_total_pyr = [corr3_off_total_pyr corr3_off_overall{i}];
    corr4_off_total_pyr = [corr4_off_total_pyr corr4_off_overall{i}];
    corr7_off_total_pyr = [corr7_off_total_pyr corr7_off_overall{i}];
    lf3_dci_total_pyr = [lf3_dci_total_pyr lf3_dci_overall{i}];
    lf4_dci_total_pyr = [lf4_dci_total_pyr lf4_dci_overall{i}];
    lf7_dci_total_pyr = [lf7_dci_total_pyr lf7_dci_overall{i}];
    wcv_dci_total_pyr = [wcv_dci_total_pyr wcv_dci_overall{i}];
    lf3_period_total_pyr = [lf3_period_total_pyr lf3_period_overall{i}];
    lf4_period_total_pyr = [lf4_period_total_pyr lf4_period_overall{i}];
    lf7_period_total_pyr = [lf7_period_total_pyr lf7_period_overall{i}];
    wcv_period_total_pyr = [wcv_period_total_pyr wcv_period_overall{i}];
    
end

%hist correlation parameters
corr_grid = [-1:0.01:1];
pyr_corr3_amp_hist = hist(corr3_amp_total_pyr,corr_grid);
pyr_corr4_amp_hist = hist(corr4_amp_total_pyr,corr_grid);
pyr_corr7_amp_hist = hist(corr7_amp_total_pyr,corr_grid);

pyr_corr3_coeff_hist = hist(corr3_coeff_total_pyr,corr_grid);
pyr_corr4_coeff_hist = hist(corr4_coeff_total_pyr,corr_grid);
pyr_corr7_coeff_hist = hist(corr7_coeff_total_pyr,corr_grid);

stell_corr3_amp_hist = hist(corr3_amp_total_stell,corr_grid);
stell_corr4_amp_hist = hist(corr4_amp_total_stell,corr_grid);
stell_corr7_amp_hist = hist(corr7_amp_total_stell,corr_grid);

stell_corr3_coeff_hist = hist(corr3_coeff_total_stell,corr_grid);
stell_corr4_coeff_hist = hist(corr4_coeff_total_stell,corr_grid);
stell_corr7_coeff_hist = hist(corr7_coeff_total_stell,corr_grid);

pyr_corr3_amp_hist = pyr_corr3_amp_hist/sum(pyr_corr3_amp_hist);
pyr_corr4_amp_hist = pyr_corr4_amp_hist/sum(pyr_corr4_amp_hist);
pyr_corr7_amp_hist = pyr_corr7_amp_hist/sum(pyr_corr7_amp_hist);

pyr_corr3_coeff_hist = pyr_corr3_coeff_hist/sum(pyr_corr3_coeff_hist);
pyr_corr4_coeff_hist = pyr_corr4_coeff_hist/sum(pyr_corr4_coeff_hist);
pyr_corr7_coeff_hist = pyr_corr7_coeff_hist/sum(pyr_corr7_coeff_hist);

stell_corr3_amp_hist = stell_corr3_amp_hist/sum(stell_corr3_amp_hist);
stell_corr4_amp_hist = stell_corr4_amp_hist/sum(stell_corr4_amp_hist);
stell_corr7_amp_hist = stell_corr7_amp_hist/sum(stell_corr7_amp_hist);

stell_corr3_coeff_hist = stell_corr3_coeff_hist/sum(stell_corr3_coeff_hist);
stell_corr4_coeff_hist = stell_corr4_coeff_hist/sum(stell_corr4_coeff_hist);
stell_corr7_coeff_hist = stell_corr7_coeff_hist/sum(stell_corr7_coeff_hist);

pyr_corr3_off_hist = hist(corr3_off_total_pyr,corr_grid);
pyr_corr4_off_hist = hist(corr4_off_total_pyr,corr_grid);
pyr_corr7_off_hist = hist(corr7_off_total_pyr,corr_grid);

pyr_corr3_off_hist = pyr_corr3_off_hist/sum(pyr_corr3_off_hist);
pyr_corr4_off_hist = pyr_corr4_off_hist/sum(pyr_corr4_off_hist);
pyr_corr7_off_hist = pyr_corr7_off_hist/sum(pyr_corr7_off_hist);

stell_corr3_off_hist = hist(corr3_off_total_stell,corr_grid);
stell_corr4_off_hist = hist(corr4_off_total_stell,corr_grid);
stell_corr7_off_hist = hist(corr7_off_total_stell,corr_grid);

stell_corr3_off_hist = stell_corr3_off_hist/sum(stell_corr3_off_hist);
stell_corr4_off_hist = stell_corr4_off_hist/sum(stell_corr4_off_hist);
stell_corr7_off_hist = stell_corr7_off_hist/sum(stell_corr7_off_hist);



%hist dci parameters
dci_grid = [-1:.01:1];

pyr_lf3_dci_hist = hist(lf3_dci_total_pyr,dci_grid);
pyr_lf4_dci_hist = hist(lf4_dci_total_pyr,dci_grid);
pyr_lf7_dci_hist = hist(lf7_dci_total_pyr,dci_grid);

pyr_lf3_dci_hist = pyr_lf3_dci_hist/sum(pyr_lf3_dci_hist);
pyr_lf4_dci_hist = pyr_lf4_dci_hist/sum(pyr_lf4_dci_hist);
pyr_lf7_dci_hist = pyr_lf7_dci_hist/sum(pyr_lf7_dci_hist);

pyr_wcv_dci_hist = hist(wcv_dci_total_pyr,dci_grid);
pyr_wcv_dci_hist = pyr_wcv_dci_hist/sum(pyr_wcv_dci_hist);

stell_lf3_dci_hist = hist(lf3_dci_total_stell,dci_grid);
stell_lf4_dci_hist = hist(lf4_dci_total_stell,dci_grid);
stell_lf7_dci_hist = hist(lf7_dci_total_stell,dci_grid);

stell_lf3_dci_hist = stell_lf3_dci_hist/sum(stell_lf3_dci_hist);
stell_lf4_dci_hist = stell_lf4_dci_hist/sum(stell_lf4_dci_hist);
stell_lf7_dci_hist = stell_lf7_dci_hist/sum(stell_lf7_dci_hist);

stell_wcv_dci_hist = hist(wcv_dci_total_stell,dci_grid);
stell_wcv_dci_hist = stell_wcv_dci_hist/sum(stell_wcv_dci_hist);

total_lf3 = [lf3_dci_total_pyr lf3_dci_total_stell];
total_lf4 = [lf4_dci_total_pyr lf4_dci_total_stell];
total_lf7 = [lf7_dci_total_pyr lf7_dci_total_stell];

total_lf3_hist = hist(total_lf3,dci_grid);
total_lf4_hist = hist(total_lf4,dci_grid);
total_lf7_hist = hist(total_lf7,dci_grid);

total_lf3_hist = total_lf3_hist/sum(total_lf3_hist);
total_lf4_hist = total_lf4_hist/sum(total_lf4_hist);
total_lf7_hist = total_lf7_hist/sum(total_lf7_hist);


%hist period params
period_grid = [0:.01:60];
pyr_lf3_period_hist = hist(lf3_period_total_pyr,period_grid);
pyr_lf3_period_hist = pyr_lf3_period_hist/sum(pyr_lf3_period_hist);
pyr_wcv_period_hist = hist(wcv_period_total_pyr,period_grid);
pyr_wcv_period_hist = pyr_wcv_period_hist/sum(pyr_wcv_period_hist);
stell_lf3_period_hist = hist(lf3_period_total_stell,period_grid);
stell_lf3_period_hist = stell_lf3_period_hist/sum(stell_lf3_period_hist);
stell_wcv_period_hist = hist(wcv_period_total_stell,period_grid);
stell_wcv_period_hist = stell_wcv_period_hist/sum(stell_wcv_period_hist);



save overall_data_jmm_multLFP_2


% figure
% subplot(3,1,1)
% plot(corr_grid,pyr_corr3_amp_hist,'linewidth',2)
% hold on
% plot(corr_grid,stell_corr3_amp_hist,'r','linewidth',2)
% xlabel('Peak Cross Correlation','FontSize',12)
% ylabel('Probability','FontSize',12)
% legend('Pyramidal','Stellate')
% title('Max XCorr LF3','FontSize',14)
% subplot(3,1,2)
% plot(corr_grid,pyr_corr4_amp_hist,'linewidth',2)
% hold on
% plot(corr_grid,stell_corr4_amp_hist,'r','linewidth',2)
% xlabel('Peak Cross Correlation','FontSize',12)
% ylabel('Probability','FontSize',12)
% legend('Pyramidal','Stellate')
% title('Max XCorr LF4','FontSize',14)
% subplot(3,1,3)
% plot(corr_grid,pyr_corr7_amp_hist,'linewidth',2)
% hold on
% plot(corr_grid,stell_corr7_amp_hist,'r','linewidth',2)
% xlabel('Peak Cross Correlation','FontSize',12)
% ylabel('Probability','FontSize',12)
% legend('Pyramidal','Stellate')
% title('Max XCorr LF7','FontSize',14)
% 
% 
% 
% figure
% subplot(3,1,1)
% plot(corr_grid,pyr_corr3_coeff_hist,'linewidth',2)
% hold on
% plot(corr_grid,stell_corr3_coeff_hist,'r','linewidth',2)
% xlabel('Correlation','FontSize',12)
% ylabel('Probability','FontSize',12)
% legend('Pyramidal','Stellate')
% title('Correlation LF3','FontSize',14)
% subplot(3,1,2)
% plot(corr_grid,pyr_corr4_coeff_hist,'linewidth',2)
% hold on
% plot(corr_grid,stell_corr4_coeff_hist,'r','linewidth',2)
% xlabel('Correlation','FontSize',12)
% ylabel('Probability','FontSize',12)
% legend('Pyramidal','Stellate')
% title('Correlation LF4','FontSize',14)
% subplot(3,1,3)
% plot(corr_grid,pyr_corr7_coeff_hist,'linewidth',2)
% hold on
% plot(corr_grid,stell_corr7_coeff_hist,'r','linewidth',2)
% xlabel('Correlation','FontSize',12)
% ylabel('Probability','FontSize',12)
% legend('Pyramidal','Stellate')
% title('Correlation LF7','FontSize',14)
% 
% figure
% subplot(3,1,1)
% plot(corr_grid,pyr_corr3_off_hist,'linewidth',2)
% hold on
% plot(corr_grid,stell_corr3_off_hist,'r','linewidth',2)
% xlabel('Peak XCorr Offset (s)','FontSize',12)
% ylabel('Probability','FontSize',12)
% legend('Pyramidal','Stellate')
% title('Peak Xcorr Offset LF3','FontSize',14)
% subplot(3,1,2)
% plot(corr_grid,pyr_corr4_off_hist,'linewidth',2)
% hold on
% plot(corr_grid,stell_corr4_off_hist,'r','linewidth',2)
% xlabel('Peak XCorr Offset (s)','FontSize',12)
% ylabel('Probability','FontSize',12)
% legend('Pyramidal','Stellate')
% title('Peak Xcorr Offset LF4','FontSize',14)
% subplot(3,1,3)
% plot(corr_grid,pyr_corr7_off_hist,'linewidth',2)
% hold on
% plot(corr_grid,stell_corr7_off_hist,'r','linewidth',2)
% xlabel('Peak XCorr Offset (s)','FontSize',12)
% ylabel('Probability','FontSize',12)
% legend('Pyramidal','Stellate')
% title('Peak Xcorr Offset LF7','FontSize',14)
% 
% figure
% subplot(3,1,1)
% plot(dci_grid,pyr_lf3_dci_hist,'linewidth',2)
% hold on
% plot(dci_grid,stell_lf3_dci_hist,'r','linewidth',2)
% legend('Pyramidal','Stellate')
% xlabel('DCI')
% ylabel('Probability','FontSize',14)
% title('LF3 DCI Distributions','FontSize',16)
% subplot(3,1,2)
% plot(dci_grid,pyr_lf4_dci_hist,'linewidth',2)
% hold on
% plot(dci_grid,stell_lf4_dci_hist,'r','linewidth',2)
% legend('Pyramidal','Stellate')
% xlabel('DCI')
% ylabel('Probability','FontSize',14)
% title('LF4 DCI Distributions','FontSize',16)
% subplot(3,1,3)
% plot(dci_grid,pyr_lf7_dci_hist,'linewidth',2)
% hold on
% plot(dci_grid,stell_lf7_dci_hist,'r','linewidth',2)
% legend('Pyramidal','Stellate')
% xlabel('DCI')
% ylabel('Probability','FontSize',14)
% title('LF7 DCI Distributions','FontSize',16)
% 
% figure
% plot(dci_grid,pyr_wcv_dci_hist,'k','linewidth',2)
% hold on
% plot(dci_grid,stell_wcv_dci_hist,'r','linewidth',2)
% plot(dci_grid,total_lf3_hist,'g','linewidth',2)
% plot(dci_grid,total_lf4_hist,'c','linewidth',2)
% plot(dci_grid,total_lf7_hist,'b','linewidth',2)
% xlabel('DCI')
% ylabel('Probability','FontSize',14)
% title('WCV DCI Distributions','FontSize',16)
% legend('Pyramidal','Stellate','LF3','LF4','LF7')
% 
% 
% figure
% plot(period_grid,fgsmooth(pyr_wcv_period_hist,3),'linewidth',2)
% hold on
% plot(period_grid,fgsmooth(stell_wcv_period_hist,3),'r','linewidth',2)
% shg
% plot(period_grid,fgsmooth(stell_lf3_period_hist,3),'k','linewidth',2)
% xlabel('Period (Hz)')
% ylabel('Probability','FontSize',14)
% title('Period Distributions','FontSize',16)
% legend('Pyramidal WCV','Stellate WCV','AVG LFP')


%check for depth dependance
load depths

figure
plot(depth(1:17),med_period(1:17),'o')
hold on
plot(depth(18:end),med_period(18:end),'ro')
legend('stellate','pyramidal')
title('Depth v period','FontSize',14)

figure
subplot(3,1,1)
plot(depth(1:17),med_corr3_amp(1:17),'o')
hold on
plot(depth(18:end),med_corr3_amp(18:end),'ro')
legend('stellate','pyramidal')
title('Depth v peak xcorr LF3','FontSize',14)
subplot(3,1,2)
plot(depth(1:17),med_corr4_amp(1:17),'o')
hold on
plot(depth(18:end),med_corr4_amp(18:end),'ro')
legend('stellate','pyramidal')
title('Depth v peak xcorr LF4','FontSize',14)
subplot(3,1,3)
plot(depth(1:17),med_corr7_amp(1:17),'o')
hold on
plot(depth(18:end),med_corr7_amp(18:end),'ro')
legend('stellate','pyramidal')
title('Depth v peak xcorr LF7','FontSize',14)

figure
subplot(3,1,1)
plot(depth(1:17),med_corr3_coef(1:17),'o')
hold on
plot(depth(18:end),med_corr3_coef(18:end),'ro')
legend('stellate','pyramidal')
title('Depth v corr coef LF3','FontSize',14)
subplot(3,1,2)
plot(depth(1:17),med_corr4_coef(1:17),'o')
hold on
plot(depth(18:end),med_corr4_coef(18:end),'ro')
legend('stellate','pyramidal')
title('Depth v corr coef LF4','FontSize',14)
subplot(3,1,3)
plot(depth(1:17),med_corr7_coef(1:17),'o')
hold on
plot(depth(18:end),med_corr7_coef(18:end),'ro')
legend('stellate','pyramidal')
title('Depth v corr coef LF7','FontSize',14)

figure
subplot(3,1,1)
plot(depth(1:17),med_corr3_off(1:17),'o')
hold on
plot(depth(18:end),med_corr3_off(18:end),'ro')
legend('stellate','pyramidal')
title('Depth v corr offset LF3','FontSize',14)
subplot(3,1,2)
plot(depth(1:17),med_corr4_off(1:17),'o')
hold on
plot(depth(18:end),med_corr4_off(18:end),'ro')
legend('stellate','pyramidal')
title('Depth v corr offset LF4','FontSize',14)
subplot(3,1,3)
plot(depth(1:17),med_corr7_off(1:17),'o')
hold on
plot(depth(18:end),med_corr7_off(18:end),'ro')
legend('stellate','pyramidal')
title('Depth v corr offset LF7','FontSize',14)

figure
plot(depth(1:17),med_wcv_dci(1:17),'o')
hold on
plot(depth(18:end),med_wcv_dci(18:end),'ro')
legend('stellate','pyramidal')
title('Depth v DCI','FontSize',14)

