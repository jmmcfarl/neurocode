
clear all

load E:\WC_Germany\JMM_Analysis\file_directories_jmm.mat
% load E:\WC_Germany\JMM_Analysis\overall_data_jmm.mat

%initialize overall variables
    lfp_dci_overall = cell(1,length(dir_array));
    wcv_dci_overall = cell(1,length(dir_array));
    diff_dci_overall = cell(1,length(dir_array));
    corr_amp_overall = cell(1,length(dir_array));
    corr_off_overall = cell(1,length(dir_array));
    corr_coeff_overall = cell(1,length(dir_array));
    lf8_period_overall = cell(1,length(dir_array));
    wcv_period_overall = cell(1,length(dir_array));
% 
    d = 30;

%cycle through directory tree
while d <= length(dir_array)
    
    cd(dir_array{d})
    pwd


    if ~exist('./used_data.mat')
        %calculate data
        hip_wc_lfp_spk_shift_jmm
    else
        load used_data
    end

    %make sure there is a folder for new figures
    if ~exist('./jmm_new')
        !mkdir jmm_new
    end


    %zscore
    wcv_z = wcv_minus_spike - mean(wcv_minus_spike);
    lf8_z = lf8 - mean(lf8);
    wcv_z = wcv_z/std(wcv_z);
    lf8_z = lf8_z/std(lf8_z);



    windowSize = 20; %window size in s
    windowSizeBins = round(windowSize/dt);
    timeSlide = 2; %amount to slide time window in s
    binSlide = round(timeSlide/dt);
    numWins = round((length(lf8)-windowSizeBins)/binSlide); %number of windows
    maxLag = round(1/dt); %max lag for cross cov calc
    Fs = mean(CSC8_SampleFrequencies);
    %initializations
    lfp_dci = zeros(1,numWins);
    wcv_dci = lfp_dci;
    dci_ratio = lfp_dci;
    corrAmp = lfp_dci;
    corrOff = lfp_dci;
    corrAbsAmp = lfp_dci;
    corrAbsOff = lfp_dci;
    corrCoSim = lfp_dci;

    for i = 1:numWins

        windowCent = (i-1)*binSlide+round(windowSizeBins/2);
        windowBeg = windowCent - round(windowSizeBins/2)+1;
        windowEnd = windowCent + round(windowSizeBins/2);
        lfpSeg = lf8_z(windowBeg:windowEnd);
        %     wcvSeg = wcv(windowBeg:windowEnd);
        wcvSeg_MS = wcv_z(windowBeg:windowEnd);

        lof = .1; hif = 2; nyqf = Fs/2; %set bandpass range
        lof = lof/nyqf; hif = hif/nyqf;
        [b,a] = butter(2, [lof hif]);

        %filter data
        lfpSeg = filtfilt(b,a,lfpSeg);
        wcvSeg_MS = filtfilt(b,a,wcvSeg_MS);

        %calculate duty cycle index
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

        temp = corrcoef(lfpSeg,wcvSeg_MS);
        if ~isempty(temp)
            corrCoSim(i) = temp(2,1);
        else
            corrCoSim(i) = NaN;
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


    lfp_dci_overall{d} = lfp_dci;
    wcv_dci_overall{d} = wcv_dci;
    diff_dci_overall{d} = dci_ratio;
    corr_amp_overall{d} = corrAmp;
    corr_off_overall{d} = corrOff;
    corr_coeff_overall{d} = corrCoSim;
    lf8_period_overall{d} = lf8_period;
    wcv_period_overall{d} = wcv_period;

   
    

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
% 
% 
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
% 
% 
%     figure
%     set(gcf,'PaperUnits','centimeters');
%     set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
%     set(gcf,'PaperPosition',[0,0,(get(gcf,'PaperSize'))])
%     subplot(2,1,1)
%     plot(timeAxis,corrAmp,'linewidth',2)
%     hold on
%     plot(timeAxis,corrCoSim,'r','linewidth',2)
%     legend('Max XCorr','Corr Coef')
%     xlabel('Time (s)','FontSize',12)
%     ylabel('Correlation Coefficient','FontSize',12)
%     title('Correlation Peak Amp vs Time','FontSize',14)
%     subplot(2,1,2)
%     plot(timeAxis,corrOff,'linewidth',2)
%     xlabel('Time (s)','FontSize',12)
%     ylabel('Peak Correlation Offset','FontSize',12)
%     title('Correlation Peak Offset vs Time','FontSize',14)
%     tname = ['./jmm_new/Correlation_xandco_' f_names{d}];
%     print('-dpng',tname)
%     tname = ['E:\WC_Germany\JMM_Analysis\corr_xandco\' f_names{d}];
%     print('-dpng',tname)
%     close all

        d = d +1
        
%             save E:\WC_Germany\JMM_Analysis\overall_data_jmm.mat *overall d

    %prepare for next iteration
    clear all

    load E:\WC_Germany\JMM_Analysis\file_directories_jmm.mat
    load E:\WC_Germany\JMM_Analysis\overall_data_jmm.mat

%     d
    
end


cd 'E:\WC_Germany\JMM_Analysis'

%calculate median values of distributions
for i = 1:43
    med_period(i) = median(wcv_period_overall{i});
    med_corr_amp(i) = median(corr_amp_overall{i});
    med_corr_coef(i) = median(corr_coeff_overall{i})
    med_corr_off(i) = median(corr_off_overall{i});
    med_wcv_dci(i) = median(wcv_dci_overall{i});
    median_diff_dci(i) = median(diff_dci_overall{i});
end

%initialize total distribution vectors
corr_amp_total_stell = [];
corr_coeff_total_stell = [];
corr_off_total_stell = [];
diff_dci_total_stell = [];
lfp_dci_total_stell = [];
wcv_dci_total_stell = [];
lf8_period_total_stell = [];
wcv_period_total_stell = [];
corr_amp_total_pyr = [];
corr_coeff_total_pyr = [];
corr_off_total_pyr = [];
diff_dci_total_pyr = [];
lfp_dci_total_pyr = [];
wcv_period_total_pyr = [];
lf8_period_total_pyr = [];
wcv_dci_total_pyr = [];

%for stellate
for i = 1:17
    corr_amp_total_stell = [corr_amp_total_stell corr_amp_overall{i}];
    corr_coeff_total_stell = [corr_coeff_total_stell corr_coeff_overall{i}];
    corr_off_total_stell = [corr_off_total_stell corr_off_overall{i}];
    diff_dci_total_stell = [diff_dci_total_stell diff_dci_overall{i}];
    lfp_dci_total_stell = [lfp_dci_total_stell lfp_dci_overall{i}];
    wcv_dci_total_stell = [wcv_dci_total_stell wcv_dci_overall{i}];
    lf8_period_total_stell = [lf8_period_total_stell lf8_period_overall{i}];
    wcv_period_total_stell = [wcv_period_total_stell wcv_period_overall{i}];
end
%for pyr
for i = 18:42
    corr_amp_total_pyr = [corr_amp_total_pyr corr_amp_overall{i}];
    corr_coeff_total_pyr = [corr_coeff_total_pyr corr_coeff_overall{i}];
    corr_off_total_pyr = [corr_off_total_pyr corr_off_overall{i}];
    diff_dci_total_pyr = [diff_dci_total_pyr diff_dci_overall{i}];
    lfp_dci_total_pyr = [lfp_dci_total_pyr lfp_dci_overall{i}];
    wcv_dci_total_pyr = [wcv_dci_total_pyr wcv_dci_overall{i}];
    lf8_period_total_pyr = [lf8_period_total_pyr lf8_period_overall{i}];
    wcv_period_total_pyr = [wcv_period_total_pyr wcv_period_overall{i}];
end

%hist correlation parameters
corr_grid = [-1:0.01:1];
pyr_corr_amp_hist = hist(corr_amp_total_pyr,corr_grid);
pyr_corr_coeff_hist = hist(corr_coeff_total_pyr,corr_grid);
stell_corr_amp_hist = hist(corr_amp_total_stell,corr_grid);
stell_corr_coeff_hist = hist(corr_coeff_total_stell,corr_grid);
pyr_corr_amp_hist = pyr_corr_amp_hist/sum(pyr_corr_amp_hist);
pyr_corr_coeff_hist = pyr_corr_coeff_hist/sum(pyr_corr_coeff_hist);
stell_corr_amp_hist = stell_corr_amp_hist/sum(stell_corr_amp_hist);
stell_corr_coeff_hist = stell_corr_coeff_hist/sum(stell_corr_coeff_hist);

%hist dci parameters
dci_grid = [-1:.01:1];
pyr_lfp_dci_hist = hist(lfp_dci_total_pyr,dci_grid);
pyr_lfp_dci_hist = pyr_lfp_dci_hist/sum(pyr_lfp_dci_hist);
pyr_wcv_dci_hist = hist(wcv_dci_total_pyr,dci_grid);
pyr_wcv_dci_hist = pyr_wcv_dci_hist/sum(pyr_wcv_dci_hist);
pyr_diff_dci_hist = hist(diff_dci_total_pyr,dci_grid);
pyr_diff_dci_hist = pyr_diff_dci_hist/sum(pyr_diff_dci_hist);
stell_lfp_dci_hist = hist(lfp_dci_total_stell,dci_grid);
stell_lfp_dci_hist = stell_lfp_dci_hist/sum(stell_lfp_dci_hist);
stell_wcv_dci_hist = hist(wcv_dci_total_stell,dci_grid);
stell_wcv_dci_hist = stell_wcv_dci_hist/sum(stell_wcv_dci_hist);
stell_diff_dci_hist = hist(diff_dci_total_stell,dci_grid);
stell_diff_dci_hist = stell_diff_dci_hist/sum(stell_diff_dci_hist);

total_lfp = [lfp_dci_total_pyr lfp_dci_total_stell];
total_lfp_hist = hist(total_lfp,dci_grid);
total_lfp_hist = total_lfp_hist/sum(total_lfp_hist);

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




plot(corr_grid,pyr_corr_amp_hist,'linewidth',2)
hold on
plot(corr_grid,stell_corr_amp_hist,'r','linewidth',2)
shg
xlabel('Peak Cross Correlation','FontSize',14)
ylabel('Probability','FontSize',14)
legend('Pyramidal','Stellate')
xlim([-0.5 1])
xlim([-0.3 1])
title('Cross Correlation Peak Distributions','FontSize',16)



plot(corr_grid,pyr_corr_coeff_hist,'linewidth',2)
hold on
plot(corr_grid,stell_corr_coeff_hist,'r','linewidth',2)
pyr_corr_coeff_hist = pyr_corr_coeff_hist/sum(pyr_corr_coeff_hist);
plot(corr_grid,pyr_corr_coeff_hist,'linewidth',2)
hold on
plot(corr_grid,stell_corr_coeff_hist,'r','linewidth',2)
xlabel('Correlation','FontSize',14)
ylabel('Probability','FontSize',14)
legend('Pyramidal','Stellate')
title('Correlation Distributions','FontSize',16)



plot(dci_grid,pyr_lfp_dci_hist,'linewidth',2)
hold on
plot(dci_grid,stell_lfp_dci_hist,'r','linewidth',2)
legend('Pyramidal','Stellate')
xlabel('DCI')
ylabel('Probability','FontSize',14)
title('LFP DCI Distributions','FontSize',16)
plot(dci_grid,pyr_wcv_dci_hist,'linewidth',2)
hold on
plot(dci_grid,stell_wcv_dci_hist,'r','linewidth',2)
legend('Pyramidal','Stellate')
xlabel('DCI')
ylabel('Probability','FontSize',14)
title('WCV DCI Distributions','FontSize',16)
plot(dci_grid,total_lfp_hist,'k','linewidth',2)
legend('Pyramidal','Stellate','AVG LFP')


plot(dci_grid,pyr_diff_dci_hist,'linewidth',2)
hold on
plot(dci_grid,stell_diff_dci_hist,'r','linewidth',2)
shg
legend('Pyramidal','Stellate')
xlabel('DCI Difference (WCV - LFP)')
ylabel('Probability','FontSize',14)
title('DCI Difference (WCV-LFP) Distributions','FontSize',16)


plot(period_grid,pyr_lf8_period_hist,'linewidth',2)
hold on
plot(period_grid,stell_lf8_period_hist,'r','linewidth',2)
shg
plot(period_grid,fgsmooth(pyr_lf8_period_hist,3),'linewidth',2)
hold on
plot(period_grid,fgsmooth(stell_lf8_period_hist,3),'r','linewidth',2)
shg
figure
plot(period_grid,fgsmooth(pyr_wcv_period_hist,3),'linewidth',2)
hold on
plot(period_grid,fgsmooth(stell_wcv_period_hist,3),'r','linewidth',2)
shg
plot(period_grid,fgsmooth(stell_lf8_period_hist,3),'k','linewidth',2)
xlabel('Period (Hz)')
ylabel('Probability','FontSize',14)
title('Period Distributions','FontSize',16)
legend('Pyramidal WCV','Stellate WCV','AVG LFP')
save overall_data_jmm