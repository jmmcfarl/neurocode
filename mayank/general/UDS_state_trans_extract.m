%% GET UP STATE AMPLITUDE BATCH

clear all

load E:\WC_Germany\JMM_Analysis\file_directories_jmm.mat


%initialize overall variables
median_period = zeros(1,43);
period_vals = cell(1,43);

d = 1;

while d <= length(dir_array)

    cd(dir_array{d})
    pwd

    load used_data CSC8_SampleFrequencies wcv_minus_spike dt synct

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
    lof = .1; hif = 20; nyqf = Fs/2;
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
    down_medwcv = downsample(medwcv,dsf);

    [peakslope_up_wcv,maxslope_up_wcv] = findpeaks(down_diff_lofwcv,'minpeakheight',0.5,'minpeakdistance',round(0.2/dt/dsf));
    sm_down_diff_medwcv = -fgsmooth(abs(down_diff_medwcv),10);

    % go through all up slope peaks
    max_dist = round(2/dt/dsf);
    upThresh = -0.2;
    wcv_up_trans_points = [];
    max_val = -0.6;
    look_ahead = 20; %number of downsampled points to look past max slope point for up state trans
    
    for i = 1:length(maxslope_up_wcv)

        %find preceding point where slope was below some threshold and
        %signal was below max val
        tempUpPoint = dsf*find((sm_down_diff_medwcv(1:(maxslope_up_wcv(i)+look_ahead))>upThresh)&(down_medwcv(1:(maxslope_up_wcv(i)+look_ahead))<max_val),1,'last');


        %%make sure that the signal is below mean at this point
%         if medwcv(tempUpPoint) > -0.5
%             tempUpPoint = [];
%         end
        %make sure that this point is within some distance of the slope peak
        %         elseif abs(up_start_temp - maxslope_up_wcv(i)) > max_dist
        %             up_start_temp = [];
        %         end

        wcv_up_trans_points = [wcv_up_trans_points tempUpPoint];

    end


%     downThresh = -0.3;
%     wcv_down_trans_points = [];
%     minDiff = round(Fs/5);
%     
%     badUpPoints = [];
%     for i = 1:length(wcv_up_trans_points)
% 
% 
%         tempDownPoint = wcv_up_trans_points(i)+minDiff+dsf*find((sm_down_diff_medwcv(round((wcv_up_trans_points(i)+minDiff)/dsf):end)>downThresh)&(down_medwcv(round((wcv_up_trans_points(i)+minDiff)/dsf):end)<max_val),1,'first');
% 
%         %          if medwcv(tempDownPoint) > -0.5
%         %              tempDownPoint = [];
%         %          end
% 
%         if isempty(tempDownPoint)
%             badUpPoints = [badUpPoints i];
%         end
%         
%         wcv_down_trans_points = [wcv_down_trans_points tempDownPoint];
% 
%     end

%     wcv_up_trans_points(badUpPoints) = [];
    
%     wcv_down_trans_points(wcv_down_trans_points > length(synct)) = [];


    %make sure first transition is up
%     if min(wcv_down_trans_points) < min(wcv_up_trans_points)
%         wcv_down_trans_points(1) = [];
%     end
%     %make sure transitions alternate
%     for i = 1:length(wcv_up_trans_points)-1
%         proceed = 0;
%         while ~proceed
%             if wcv_up_trans_points(i+1)<wcv_down_trans_points(i)
%                 wcv_up_trans_points(i+1) = [];
%             else 
%                 proceed = 1;
%             end
%         end
%         proceed = 0;
%         while ~proceed
%             if wcv_down_trans_points(i+1)< wcv_up_trans_points(i+1)
%                 wcv_down_trans_points(i+1) = [];
%             else
%                 proceed = 1;
%             end
%         end
%     end
    %make sure end on down transition
%     if max(wcv_up_trans_points) > max(wcv_down_trans_points)
%         wcv_up_trans_points(end) = [];
%     end
    
%     %create vector to indicate state
%     isUp = zeros(size(wcv_minus_spike));
%     if length(wcv_up_trans_points) ~= length(wcv_down_trans_points)
%         disp('ERROR DETECTING TRANSITIONS')
%     end
%     for i = 1:length(wcv_up_trans_points)
%        
%         isUp(wcv_up_trans_points(i):wcv_down_trans_points(i)) = true;
%         
%     end


    %calculate DCI distribution in segmented windows
    windowSize = 20; %window size in s
    windowSizeBins = round(windowSize/dt);
    timeSlide = 2; %amount to slide time window in s
    binSlide = round(timeSlide/dt);
    numWins = round((length(wcv_minus_spike)-windowSizeBins)/binSlide); %number of windows

    cur_period = zeros(1,numWins);
    for i = 1:numWins

        windowCent = (i-1)*binSlide+round(windowSizeBins/2);
        windowBeg = windowCent - round(windowSizeBins/2)+1;
        windowEnd = windowCent + round(windowSizeBins/2);

        cur_period(i) = length(find((wcv_up_trans_points > windowBeg) & wcv_up_trans_points < windowEnd))/windowSize;

    end
    
    median_period(d) = median(cur_period);
    period_vals{d} = cur_period;
    
    d = d +1

    save E:\WC_Germany\JMM_Analysis\overall_data_jmm_wcv_UDS_trans_period.mat period_vals d median*

    %prepare for next iteration
    clear all

    load E:\WC_Germany\JMM_Analysis\file_directories_jmm.mat
    load E:\WC_Germany\JMM_Analysis\overall_data_jmm_wcv_UDS_trans_period.mat period_vals d median*


end




% dci_grid = [-1:.01:1];
% dci_mat = zeros(43,length(dci_grid));
% for i = 1:43
%     dci_mat(i,:) = hist(duty_cycle_vals{i},dci_grid);
% end
% 
% figure
% pcolor(dci_grid,[1:43],dci_mat);shading flat;colorbar
% xlabel('DCI','FontSize',14)
% ylabel('Cell Index','FontSize',14)
% 
% load depths
% figure
% plot(depth(1:17),median_duty_cycle(1:17),'o')
% hold on
% plot(depth(18:end),median_duty_cycle(18:end),'ro')
% legend('stellate','pyramidal')
% 
% figure
% plot(depth(1:17),overall_duty_cycle(1:17),'o')
% hold on
% plot(depth(18:end),overall_duty_cycle(18:end),'ro')
% legend('stellate','pyramidal')
% 
% 
% 
% 
% 
% 
