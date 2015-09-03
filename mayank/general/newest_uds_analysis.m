%% new wcv analysis
%cell array of directories to scan
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


dci_grid = [-1:0.01:1];
lfp_dci_overall = [];
wcv_dci_overall = [];
diff_dci_overall = [];
corr_amp_overall = [];
corr_off_overall = [];
corr_abs_amp_overall = [];
corr_abs_off_overall = [];

for d = 8:13

    %change into target directory
    cd(dir_array{d})
    pwd

    %make sure there is a folder for new figures
    if ~exist('./jmm_new')
        !mkdir jmm_new
    end

    
    %load data file
    load u
    
    %take 20 second segments

    windowSize = 20; %window size in s
    windowSizeBins = round(windowSize/dt);
    timeSlide = 2; %amount to slide time window in s
    binSlide = round(timeSlide/dt);
    numWins = round((length(loflf8)-windowSizeBins)/binSlide); %number of windows
    maxLag = round(1/dt); %max lag for cross cov calc
%     NFFT = 2^13;
    Fs = mean(CSC8_SampleFrequencies);
    %initializations
    lfp_dci = zeros(1,numWins);
    wcv_dci = lfp_dci;
    dci_ratio = lfp_dci;
    corrAmp = lfp_dci;
    corrOff = lfp_dci;
    corrAbsAmp = lfp_dci;
    corrAbsOff = lfp_dci;

    for i = 1:numWins
        
        windowCent = (i-1)*binSlide+round(windowSizeBins/2);
        windowBeg = windowCent - round(windowSizeBins/2)+1;
        windowEnd = windowCent + round(windowSizeBins/2);
        lfpSeg = lf8(windowBeg:windowEnd);
        %     wcvSeg = wcv(windowBeg:windowEnd);
        wcvSeg_MS = wcv_minus_spike(windowBeg:windowEnd);

        lof = .1; hif = 2; nyqf = Fs/2; %set bandpass range
        lof = lof/nyqf; hif = hif/nyqf;
        [b,a] = butter(2, [lof hif]);

        %filter data
        lfpSeg = filtfilt(b,a,lfpSeg);
        wcvSeg_MS = filtfilt(b,a,wcvSeg_MS);

        %detrend data
%         lfpSeg = detrend(lfpSeg);
%         wcvSeg_MS = detrend(wcvSeg_MS);

        %z-score data
        lfpSeg = lfpSeg - mean(lfpSeg);
        lfpSeg = lfpSeg/std(lfpSeg);
        wcvSeg_MS = wcvSeg_MS - mean(wcvSeg_MS);
        wcvSeg_MS = wcvSeg_MS/std(wcvSeg_MS);

        %calculate duty cycle index

%         lfp_up_time = length(find(lfpSeg>0.5));
%         lfp_down_time = length(find(lfpSeg<-0.5));
%         wcv_up_time = length(find(wcvSeg_MS>0.5));
%         wcv_down_time = length(find(wcvSeg_MS<-0.5));
        lfp_up_time = length(find(lfpSeg>0));
        lfp_down_time = length(find(lfpSeg<-0));
        wcv_up_time = length(find(wcvSeg_MS>0));
        wcv_down_time = length(find(wcvSeg_MS<-0));

        lfp_dci(i) = (lfp_up_time-lfp_down_time)/(lfp_up_time+lfp_down_time);
        wcv_dci(i) = (wcv_up_time-wcv_down_time)/(wcv_up_time+wcv_down_time);

        dci_ratio(i) = wcv_dci(i) - lfp_dci(i);

        %calculate covariance peak and offset

        [corrVec,corrT] = xcov(lfpSeg,wcvSeg_MS,maxLag,'coef');
        [peakAmp,peakLoc] = max(corrVec);
        [peakAbsAmp,peakAbsLoc] = max(abs(corrVec));
        corrAmp(i) = peakAmp;
        corrOff(i) = (peakLoc - maxLag)*dt;
        corrAbsAmp(i) = peakAbsAmp;
        corrAbsOff(i) = (peakAbsLoc - maxLag)*dt;
        %

        %calculate power spectrum over data segs
%         [lfp_Pxx,F] = pwelch(lfpSeg,[],[],NFFT,Fs);
%         [wcv_Pxx,F] = pwelch(wcvSeg_MS,[],[],NFFT,Fs);

        %disregard power over 100 Hz
%         cutoff = find(F > 100,1,'first');
%         lfp_Pow(i,:) = lfp_Pxx(1:cutoff);
%         wcv_Pow(i,:) = wcv_Pxx(1:cutoff);
%         rel_Pow(i,:) = wcv_Pow(i,:)./lfp_Pow(i,:);
% 
        sprintf('%d of %d windows',i,numWins)

    end

    
            lfp_dci_overall = [lfp_dci_overall lfp_dci];
        wcv_dci_overall = [wcv_dci_overall wcv_dci];
        diff_dci_overall = [diff_dci_overall dci_ratio];

        corr_amp_overall = [corr_amp_overall corrAmp];
        corr_off_overall = [corr_off_overall corrOff];
        corr_abs_amp_overall = [corr_abs_amp_overall corrAbsAmp];
        corr_abs_off_overall = [corr_abs_off_overall corrAbsOff];
        
    %get rid of extra frequency points
%     F(cutoff+1:end) = [];

    timeAxis = [1:numWins]*timeSlide;

    % % to better view correlation offset changes set outliers as NaNs
    % zcorrOff = corrOff - mean(corrOff);
    % zcorrOff = zcorrOff/std(corrOff);
    % corrOff(abs(zcorrOff)>3) = NaN;

%     figure
%     set(gcf,'PaperUnits','centimeters');
%     set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
%     set(gcf,'PaperPosition',[0,0,(get(gcf,'PaperSize'))])
%     plot(timeAxis,lfp_dci,'linewidth',2)
%     hold on
%     plot(timeAxis,wcv_dci,'r','linewidth',2)
%     plot(timeAxis,dci_ratio,'k','linewidth',2)
%     legend('LFP DCI','WCV DCI','DCI Ratio Index')
%     xlabel('Time (s)','FontSize',14);
%     ylabel('Duty Cycle Index','FontSize',14)
%     title('Duty Cycle Index Vs Time ','FontSize',16)
%     tname = ['./jmm_new/DCI_ratio_sub_' f_names{d}];
%     print('-dpng',tname)
%     tname = ['E:\WC_Germany\JMM_Analysis\duty_cycle_index_sub\' f_names{d}];
%     print('-dpng',tname)
%     close all
clear lfp_dci wcv_dci dci_ratio
    
    
%     figure
%     set(gcf,'PaperUnits','centimeters');
%     set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
%     set(gcf,'PaperPosition',[0,0,(get(gcf,'PaperSize'))])
%     subplot(2,1,1)
%     plot(timeAxis,corrAmp,'linewidth',2)
%     hold on
%     plot(timeAxis,corrAbsAmp,'r','linewidth',2)
%     legend('Corr','Abs Corr')
%     xlabel('Time (s)','FontSize',12)
%     ylabel('Correlation Coefficient','FontSize',12)
%     title('Correlation Peak Amp vs Time','FontSize',14)
%     subplot(2,1,2)
%     plot(timeAxis,corrOff,'linewidth',2)
%     hold on
%     plot(timeAxis,corrAbsOff,'r','linewidth',2)
%     legend('Peak Corr Off','Peak Abs Corr Off')
%     xlabel('Time (s)','FontSize',12)
%     ylabel('Peak Correlation Offset','FontSize',12)
%     title('Correlation Peak Offset vs Time','FontSize',14)
%     tname = ['./jmm_new/Correlation_amp_off_abs_' f_names{d}];
%     print('-dpng',tname)
%     tname = ['E:\WC_Germany\JMM_Analysis\corr_amp_off_abs\' f_names{d}];
%     print('-dpng',tname)
%     close all
clear corrAmp corrOff
% 
%     figure
%     set(gcf,'PaperUnits','centimeters');
%     set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
%     set(gcf,'PaperPosition',[0,0,(get(gcf,'PaperSize'))])
%     subplot(2,1,1)
%     pcolor(timeAxis,F,log(lfp_Pow'));shading flat;
%     caxis([-20 0]);colorbar
%     xlabel('Time (s)','FontSize',12)
%     ylabel('Frequency (Hz)','FontSize',12)
%     title('Lfp Power','FontSize',14)
%     subplot(2,1,2)
%     pcolor(timeAxis,F,log(wcv_Pow'));shading flat;
%     caxis([-20 0]);colorbar
%     xlabel('Time (s)','FontSize',12)
%     ylabel('Frequency (Hz)','FontSize',12)
%     title('Wcv Power','FontSize',14)
%     tname = ['./jmm_new/specgram_lfp_wcv_' f_names{d}];
%     print('-dpng',tname)
%     tname = ['E:\WC_Germany\JMM_Analysis\power_spectra\' f_names{d}];
%     print('-dpng',tname)
%     close all
% clear lfp_Pow wcv_Pow
%     
%     figure
%     set(gcf,'PaperUnits','centimeters');
%     set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
%     set(gcf,'PaperPosition',[0,0,(get(gcf,'PaperSize'))])
%     pcolor(timeAxis,F,log(rel_Pow'));shading flat;colorbar
%     xlabel('Time (s)','FontSize',12)
%     ylabel('Frequency (Hz)','FontSize',12)
%     title('Relative Power (Wcv/Lfp)','FontSize',14)
%     tname = ['./jmm_new/specgram_rel_' f_names{d}];
%     print('-dpng',tname)
%     tname = ['E:\WC_Germany\JMM_Analysis\rel_power_spectra\' f_names{d}];
%     print('-dpng',tname)
%     close all
% clear rel_Pow timeAxis 

end

% pyr_corr_abs_amp_hist = corr_abs_amp_hist/sum(corr_abs_amp_hist);
% pyr_corr_abs_off_hist = corr_abs_off_hist/sum(corr_abs_off_hist);
% pyr_corr_amp_hist = corr_amp_hist/sum(corr_amp_hist);
% pyr_corr_off_hist = corr_off_hist/sum(corr_off_hist);
% pyr_lfp_dci_hist = lfp_dci_hist/sum(lfp_dci_hist);
% pyr_wcv_dci_hist = wcv_dci_hist/sum(wcv_dci_hist);
% pyr_rat_dci_hist = rat_dci_hist/sum(rat_dci_hist);

cd 'E:\WC_Germany\JMM_Analysis'
save stell_hist_data_new *overall