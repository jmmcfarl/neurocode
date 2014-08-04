%% GET UP STATE AMPLITUDE BATCH

clear all

load E:\WC_Germany\JMM_Analysis\file_directories_jmm.mat


%initialize overall variables
overall_UpStateAmp = cell(1,43);

d = 1;

while d <= length(dir_array)
    
    cd(dir_array{d})
    pwd
    
    load used_data

    Fs = mean(CSC8_SampleFrequencies);
    
    %zscore 
    wcv_z = wcv_minus_spike - mean(wcv_minus_spike);
    wcv_z = wcv_z/std(wcv_z);
    
    %filter data
    lof = .1; hif = 0.5; nyqf = Fs/2;
    lof = lof/nyqf; hif = hif/nyqf;
    [b,a] = butter(2, [lof hif]);
    lofwcv = filtfilt(b,a,wcv_minus_spike);

    %higher range bandpass
    lof = .1; hif = 10; nyqf = Fs/2;
    lof = lof/nyqf; hif = hif/nyqf;
    [b,a] = butter(2, [lof hif]);
    medwcv = filtfilt(b,a,wcv_z);


    %take derivative of filtered signals
    diff_lofwcv = [0;diff(lofwcv)];
    diff_medwcv = [0;diff(medwcv)];

    %zscore derivative signal
    diff_lofwcv = diff_lofwcv - mean(diff_lofwcv);
    diff_lofwcv = diff_lofwcv/std(diff_lofwcv);

    diff_medwcv = diff_medwcv - mean(diff_medwcv);
    diff_medwcv = diff_medwcv/std(diff_medwcv);


    %downsample and find peaks
    dsf = 20;
    down_diff_lofwcv = downsample(diff_lofwcv,dsf);
    down_diff_medwcv = downsample(diff_medwcv,dsf);

    [peakslope_up_wcv,maxslope_up_wcv] = findpeaks(down_diff_lofwcv,'minpeakheight',1,'minpeakdistance',round(0.2/dt/dsf));
    
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

    
    upStateAmp = zeros(1,length(wcv_up_trans_points));
    for i = 1:length(wcv_up_trans_points)
        
       nextUpStateAmp = lofwcv(wcv_up_trans_points(i)+find(diff_lofwcv(wcv_up_trans_points(i):end)<0,1,'first'));
       nextDownStateAmp = lofwcv(find(diff_lofwcv(1:wcv_up_trans_points(i))<0,1,'last')); 
       
       if isempty(nextUpStateAmp) | isempty(nextDownStateAmp)
           upStateAmp(i) = NaN;
       else     
       upStateAmp(i) = nextUpStateAmp - nextDownStateAmp;
       end
       
    end
    
    overall_UpStateAmp{d} = upStateAmp;
    
    d = d +1

    save E:\WC_Germany\JMM_Analysis\overall_data_jmm_wcv_amp.mat overall* d

    %prepare for next iteration
    clear all

    load E:\WC_Germany\JMM_Analysis\file_directories_jmm.mat
    load E:\WC_Germany\JMM_Analysis\overall_data_jmm_wcv_amp.mat

    
end


%% Extract median amplitude for each cell

wcv_Gains = [100 50 100 100 100 150 100 150 120 150 120 100 100 100 100 110 110 ... 
    100 50 50 40 100 100 80 80 90 100 100 90 120 100 80 70 80 80 100 80 80 100 90 90 100 100];

load depths
stellate_amp = zeros(1,17);
pyr_amp = zeros(1,26);

for i = 1:17
    
    stellate_amp(i) = nanmedian(overall_UpStateAmp{i});
    
end

for i = 18:43
    
    pyr_amp(i-17) = nanmedian(overall_UpStateAmp{i});
    
end

stellate_amp = stellate_amp./wcv_Gains(1:17);
pyr_amp = pyr_amp./wcv_Gains(18:end);

% average data points that are the same cell, diff sessions
pairIDs = [15,17,23,25,35,38,41,43];
stellate_amp(14) = mean([stellate_amp(14) stellate_amp(15)]);
stellate_amp(16) = mean([stellate_amp(16) stellate_amp(17)]);
pyr_amp(22) = mean([pyr_amp(22) pyr_amp(23)]);
pyr_amp(24) = mean([pyr_amp(24) pyr_amp(25)]);
pyr_amp(34) = mean([pyr_amp(34) pyr_amp(35)]);
pyr_amp(37) = mean([pyr_amp(37) pyr_amp(38)]);
pyr_amp(40) = mean([pyr_amp(40) pyr_amp(41)]);
pyr_amp(42) = mean([pyr_amp(42) pyr_amp(43)]);


depth(pairIDs) = [];
stellate_amp([15 17]) = [];
pyr_amp(pairIDs(3:end)-17) = [];


% load DCV_Axis

figure
plot(depth(1:15),stellate_amp,'o','MarkerSize',6,'linewidth',2)
hold on
plot(depth(16:end),pyr_amp,'ro','MarkerSize',6,'linewidth',2)
xlabel('Depth (um)','FontSize',14)
ylabel('Up State Amplitude (mV)','FontSize',14)
title('Up State Amp V Depth','FontSize',16)




    