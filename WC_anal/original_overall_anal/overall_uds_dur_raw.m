% clear all

load('C:\WC_Germany\overall_calcs\overall_dir.mat')
load C:\WC_Germany\overall_calcs\run_hist_thresh\thresh_data 

segDur = 30;
winSlide = 1;
dsf = 8;
Fsd = 2016/dsf;

thresh_smooth = 10;

niqf = 1/2;
[b,a] = butter(2,0.01/niqf,'low');
for d = 1:length(over_dir)

    cd(over_dir{d})
    pwd

    load used_data lf8 wcv_minus_spike wcv
%     load spike_time
%     if exist('spike_time.mat')
%         load spike_time
%     else
%         load spike_time_br
%     end
    
%     sm_threshold_w{d} = jmm_smooth_1d_cor(threshold_w{d},thresh_smooth);
%     sm_threshold_8{d} = jmm_smooth_1d_cor(threshold_8{d},thresh_smooth);

    sm_threshold_w{d} = filtfilt(b,a,threshold_w{d});
    sm_threshold_8{d} = filtfilt(b,a,threshold_8{d});
    
    disp('WCV')
    [up_state_dur{d},down_state_dur{d},up_trans{d},down_trans{d},sig_w{d},int_thresh_w{d}] = UDS_extract_run_hist(wcv_minus_spike,sm_threshold_w{d});
    disp('LFP')
    [up_state_dur8{d},down_state_dur8{d},up_trans8{d},down_trans8{d},sig_8{d},int_thresh_8{d}] = UDS_extract_run_hist(lf8,sm_threshold_8{d});

    % dataViewer_run_hist(sig_w{d},sig_8{d},int_thresh_w{d},int_thresh_8{d},up_trans{d},up_trans8{d},down_trans{d},down_trans8{d})


    %calculate average and max up state duration in segmented data
    nsamples = length(sig_w{d});
    numWins = floor((nsamples - segDur*Fsd)/(winSlide*Fsd));
    for w = 1:numWins

        begSamp = round(1+(w-1)*winSlide*Fsd);
        endSamp = begSamp + round(segDur*Fsd);

        %find current up states
        curUps_w = find(up_trans{d}>=begSamp&up_trans{d}<endSamp);
        curUps_8 = find(up_trans8{d}>=begSamp&up_trans8{d}<endSamp);

        curDowns_w = find(down_trans{d}>=begSamp&down_trans{d}<endSamp);
        curDowns_8 = find(down_trans8{d}>=begSamp&down_trans8{d}<endSamp);

        avg_up_dur{d}(w) = mean(up_state_dur{d}(curUps_w));
        avg_up_dur8{d}(w) = mean(up_state_dur8{d}(curUps_8));

        if ~isempty(curUps_w)
            max_up_dur{d}(w) = max(up_state_dur{d}(curUps_w));
        else
            max_up_dur{d}(w) = nan;
        end

        if ~isempty(curUps_8)
            max_up_dur8{d}(w) = max(up_state_dur8{d}(curUps_8));
        else
            max_up_dur8{d}(w) = nan;
        end
        curDowns_w(curDowns_w==length(down_trans{d})) = [];
        curDowns_8(curDowns_8==length(down_trans8{d})) = [];
        if ~isempty(curDowns_w)
            max_down_dur{d}(w) = max(down_state_dur{d}(curDowns_w));
        else
            max_down_dur{d}(w) = nan;
        end

        if ~isempty(curDowns_8)
            max_down_dur8{d}(w) = max(down_state_dur8{d}(curDowns_8));
        else
            max_down_dur8{d}(w) = nan;
        end


        avg_down_dur{d}(w) = mean(down_state_dur{d}(curDowns_w));
        avg_down_dur8{d}(w) = mean(down_state_dur8{d}(curDowns_8));

    end
    %

    des_times = find(up_state_dur8{d} > 5);
    
    end_time = length(max_up_dur{d});

    plot(max_up_dur{d},'linewidth',2)
    hold on
    plot(max_up_dur8{d},'r','linewidth',2)
    xlim([0 end_time])
    plot(up_trans8{d}(des_times)/Fsd,ones(size(des_times))*5,'k*')
    t_names = ['C:\WC_Germany\overall_calcs\UDS_dur_raw\new_method_' num2str(cell_type(d)) '_' over_names{d}]
    print('-dpng',t_names)
    close all


    save C:\WC_Germany\overall_calcs\UDS_dur_raw\UDS_raw_data_new_method up_state_dur* down_state_dur* up_trans* down_trans* 
    %
    clear lf8 wcv_minus_spike 
    
end


