clear all
cd C:\WC_Germany\DUAL-MP_recordings
load dual_mp_dir
load run_hist_thresh\UDS_dur_data

segDur = 20;
winSlide = 1;
dsf = 8;
Fsd = 2016/dsf;

thresh_smooth = 10;

niqf = 1/2;
[b,a] = butter(2,.005/niqf,'low');

for d = 1:28
    cd(dir_array{d})
    pwd

    load used_data lf*
    
    sm_threshold_8{d} = filtfilt(b,a,threshold_8{d});
     sm_threshold_13{d} = filtfilt(b,a,threshold_13{d});
    sm_threshold_15{d} = filtfilt(b,a,threshold_15{d});
    sm_threshold_16{d} = filtfilt(b,a,threshold_16{d});
   if exist('lf14')
           sm_threshold_14{d} = filtfilt(b,a,threshold_14{d});
   end
   
    disp('LF8')
    [up_state_dur8{d},down_state_dur8{d},up_trans8{d},down_trans8{d},sig_8{d},int_thresh_8{d},tot_chopped_up8(d),tot_chopped_down8(d)] = UDS_extract_run_hist(lf8,sm_threshold_8{d});
    disp('LF13')
    [up_state_dur13{d},down_state_dur13{d},up_trans13{d},down_trans13{d},sig_13{d},int_thresh_13{d},tot_chopped_up13(d),tot_chopped_down13(d)] = UDS_extract_run_hist(lf13,sm_threshold_13{d});
    disp('LF15')
    [up_state_dur15{d},down_state_dur15{d},up_trans15{d},down_trans15{d},sig_15{d},int_thresh_15{d},tot_chopped_up15(d),tot_chopped_down15(d)] = UDS_extract_run_hist(lf15,sm_threshold_15{d});
    disp('LF16')
    [up_state_dur16{d},down_state_dur16{d},up_trans16{d},down_trans16{d},sig_16{d},int_thresh_16{d},tot_chopped_up16(d),tot_chopped_down16(d)] = UDS_extract_run_hist(lf16,sm_threshold_16{d});
if exist('lf14')
        disp('LF14')
    [up_state_dur14{d},down_state_dur14{d},up_trans14{d},down_trans14{d},sig_14{d},int_thresh_14{d},tot_chopped_up14(d),tot_chopped_down14(d)] = UDS_extract_run_hist(lf14,sm_threshold_14{d});
end
%     dataViewer_run_hist(sig_w{d},sig_8{d},int_thresh_w{d},int_thresh_8{d},up_trans{d},up_trans8{d},down_trans{d},down_trans8{d})


    %calculate average and max up state duration in segmented data
    nsamples = length(sig_8{d});
    numWins = floor((nsamples - segDur*Fsd)/(winSlide*Fsd));
    for w = 1:numWins

        begSamp = round(1+(w-1)*winSlide*Fsd);
        endSamp = begSamp + round(segDur*Fsd);

        %find current up states
        curUps_8 = find(up_trans8{d}>=begSamp&up_trans8{d}<endSamp);
        curUps_13 = find(up_trans13{d}>=begSamp&up_trans13{d}<endSamp);
        curUps_15 = find(up_trans15{d}>=begSamp&up_trans15{d}<endSamp);
        curUps_16 = find(up_trans16{d}>=begSamp&up_trans16{d}<endSamp);
        if exist('lf14')
        curUps_14 = find(up_trans14{d}>=begSamp&up_trans14{d}<endSamp);
        end
        
        curDowns_8 = find(down_trans8{d}>=begSamp&down_trans8{d}<endSamp);
        curDowns_13 = find(down_trans13{d}>=begSamp&down_trans13{d}<endSamp);
        curDowns_15 = find(down_trans15{d}>=begSamp&down_trans15{d}<endSamp);
        curDowns_16 = find(down_trans16{d}>=begSamp&down_trans16{d}<endSamp);
if exist('lf14')
            curDowns_14 = find(down_trans14{d}>=begSamp&down_trans14{d}<endSamp);
end

avg_up_dur8{d}(w) = mean(up_state_dur8{d}(curUps_8));
avg_up_dur13{d}(w) = mean(up_state_dur13{d}(curUps_13));
avg_up_dur15{d}(w) = mean(up_state_dur15{d}(curUps_15));
avg_up_dur16{d}(w) = mean(up_state_dur16{d}(curUps_16));
if exist('lf14')
    avg_up_dur14{d}(w) = mean(up_state_dur14{d}(curUps_14));
end

        if ~isempty(curUps_8)
            max_up_dur8{d}(w) = max(up_state_dur8{d}(curUps_8));
        else
            max_up_dur8{d}(w) = nan;
        end
        curDowns_8(curDowns_8==length(down_trans8{d})) = [];
        if ~isempty(curDowns_8)
            max_down_dur8{d}(w) = max(down_state_dur8{d}(curDowns_8));
        else
            max_down_dur8{d}(w) = nan;
        end
        avg_down_dur8{d}(w) = mean(down_state_dur8{d}(curDowns_8));

                if ~isempty(curUps_13)
            max_up_dur13{d}(w) = max(up_state_dur13{d}(curUps_13));
        else
            max_up_dur13{d}(w) = nan;
        end
        curDowns_13(curDowns_13==length(down_trans13{d})) = [];
        if ~isempty(curDowns_13)
            max_down_dur13{d}(w) = max(down_state_dur13{d}(curDowns_13));
        else
            max_down_dur13{d}(w) = nan;
        end
        avg_down_dur13{d}(w) = mean(down_state_dur13{d}(curDowns_13));

                if ~isempty(curUps_15)
            max_up_dur15{d}(w) = max(up_state_dur15{d}(curUps_15));
        else
            max_up_dur15{d}(w) = nan;
        end
        curDowns_15(curDowns_15==length(down_trans15{d})) = [];
        if ~isempty(curDowns_15)
            max_down_dur15{d}(w) = max(down_state_dur15{d}(curDowns_15));
        else
            max_down_dur15{d}(w) = nan;
        end
        avg_down_dur15{d}(w) = mean(down_state_dur15{d}(curDowns_15));

                if ~isempty(curUps_16)
            max_up_dur16{d}(w) = max(up_state_dur16{d}(curUps_16));
        else
            max_up_dur16{d}(w) = nan;
        end
        curDowns_16(curDowns_16==length(down_trans16{d})) = [];
        if ~isempty(curDowns_16)
            max_down_dur16{d}(w) = max(down_state_dur16{d}(curDowns_16));
        else
            max_down_dur16{d}(w) = nan;
        end
        avg_down_dur16{d}(w) = mean(down_state_dur16{d}(curDowns_16));

if exist('lf14')
                   if ~isempty(curUps_14)
            max_up_dur14{d}(w) = max(up_state_dur14{d}(curUps_14));
        else
            max_up_dur14{d}(w) = nan;
        end
        curDowns_14(curDowns_14==length(down_trans14{d})) = [];
        if ~isempty(curDowns_14)
            max_down_dur14{d}(w) = max(down_state_dur14{d}(curDowns_14));
        else
            max_down_dur14{d}(w) = nan;
        end
        avg_down_dur14{d}(w) = mean(down_state_dur14{d}(curDowns_14));

 
end
    end
    %

    
    end_time = length(max_up_dur8{d});

    plot(max_up_dur8{d},'r','linewidth',2)
    hold on
    plot(max_up_dur13{d},'k','linewidth',2)
    xlim([0 end_time])
    tname = ['C:\WC_Germany\DUAL-MP_recordings\UDS_dur_raw\max_up_dur13_' f_names{d}];
    print('-dpng',tname)
    close all

    plot(max_up_dur8{d},'r','linewidth',2)
    hold on
    plot(max_up_dur15{d},'k','linewidth',2)
    xlim([0 end_time])
    tname = ['C:\WC_Germany\DUAL-MP_recordings\UDS_dur_raw\max_up_dur15_' f_names{d}];
    print('-dpng',tname)
    close all

        plot(max_up_dur8{d},'r','linewidth',2)
    hold on
    plot(max_up_dur16{d},'k','linewidth',2)
    xlim([0 end_time])
    tname = ['C:\WC_Germany\DUAL-MP_recordings\UDS_dur_raw\max_up_dur16_' f_names{d}];
    print('-dpng',tname)
    close all

    if exist('lf14')
           plot(max_up_dur8{d},'r','linewidth',2)
    hold on
    plot(max_up_dur14{d},'k','linewidth',2)
    xlim([0 end_time])
    tname = ['C:\WC_Germany\DUAL-MP_recordings\UDS_dur_raw\max_up_dur14_' f_names{d}];
    print('-dpng',tname)
    close all
 
    end
    %
    clear lf*
    
end

save C:\WC_Germany\DUAL-MP_recordings\UDS_dur_raw\UDS_raw_data up_state_dur* down_state_dur* up_trans* down_trans* 

