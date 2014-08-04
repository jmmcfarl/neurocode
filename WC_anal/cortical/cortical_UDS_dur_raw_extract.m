clear all
close all
cd C:\WC_Germany\Cortical_analysis

load C:\WC_Germany\Cortical_analysis\cortical_dir
load C:\WC_Germany\Cortical_analysis\run_hist_thresh\UDS_dur_data_w_8

segDur = 20;
winSlide = 1;
dsf = 8;
Fsd = 2016/dsf;

thresh_smooth = 10;

niqf = 1/2;
[b,a] = butter(2,.005/niqf,'low');

for d = 1:length(sess_data)

    cd(sess_data(d).directory)
    pwd

    load used_data lf8 wcv_minus_spike wcv
    
    sm_threshold_w{d} = filtfilt(b,a,threshold_w{d});
    sm_threshold_8{d} = filtfilt(b,a,threshold_8{d});
%     sm_threshold_5{d} = filtfilt(b,a,threshold_5{d});
%     sm_threshold_7{d} = filtfilt(b,a,threshold_7{d});
   
    disp('WCV')
    [up_state_dur{d},down_state_dur{d},up_trans{d},down_trans{d},sig_w{d},int_thresh_w{d},tot_chopped_up(d),tot_chopped_down(d)] = UDS_extract_run_hist(wcv_minus_spike,sm_threshold_w{d});
    disp('LF8')
    [up_state_dur8{d},down_state_dur8{d},up_trans8{d},down_trans8{d},sig_8{d},int_thresh_8{d},tot_chopped_up8(d),tot_chopped_down8(d)] = UDS_extract_run_hist(lf8,sm_threshold_8{d});
%     disp('LF5')
%     [up_state_dur5{d},down_state_dur5{d},up_trans5{d},down_trans5{d},sig_5{d},int_thresh_5{d},tot_chopped_up5(d),tot_chopped_down5(d)] = UDS_extract_run_hist(lf5,sm_threshold_5{d});
%     disp('LF7')
%     [up_state_dur7{d},down_state_dur7{d},up_trans7{d},down_trans7{d},sig_7{d},int_thresh_7{d},tot_chopped_up7(d),tot_chopped_down7(d)] = UDS_extract_run_hist(lf7,sm_threshold_7{d});
   
%     dataViewer_run_hist(sig_w{d},sig_8{d},int_thresh_w{d},int_thresh_8{d},up_trans{d},up_trans8{d},down_trans{d},down_trans8{d})


    %calculate average and max up state duration in segmented data
    nsamples = length(sig_w{d});
    numWins = floor((nsamples - segDur*Fsd)/(winSlide*Fsd));
    for w = 1:numWins

        begSamp = round(1+(w-1)*winSlide*Fsd);
        endSamp = begSamp + round(segDur*Fsd);

        %find current up states
        curUps_w = find(up_trans{d}>=begSamp&up_trans{d}<endSamp);
        curUps_8 = find(up_trans8{d}>=begSamp&up_trans8{d}<endSamp);
%         curUps_5 = find(up_trans5{d}>=begSamp&up_trans5{d}<endSamp);
%         curUps_7 = find(up_trans7{d}>=begSamp&up_trans7{d}<endSamp);

        curDowns_w = find(down_trans{d}>=begSamp&down_trans{d}<endSamp);
        curDowns_8 = find(down_trans8{d}>=begSamp&down_trans8{d}<endSamp);
%         curDowns_5 = find(down_trans5{d}>=begSamp&down_trans5{d}<endSamp);
%         curDowns_7 = find(down_trans7{d}>=begSamp&down_trans7{d}<endSamp);
        
        avg_up_dur{d}(w) = mean(up_state_dur{d}(curUps_w));
        avg_up_dur8{d}(w) = mean(up_state_dur8{d}(curUps_8));
%         avg_up_dur5{d}(w) = mean(up_state_dur5{d}(curUps_5));
%         avg_up_dur7{d}(w) = mean(up_state_dur7{d}(curUps_7));

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
%         if ~isempty(curUps_5)
%             max_up_dur5{d}(w) = max(up_state_dur5{d}(curUps_5));
%         else
%             max_up_dur5{d}(w) = nan;
%         end
%         if ~isempty(curUps_7)
%             max_up_dur7{d}(w) = max(up_state_dur7{d}(curUps_7));
%         else
%             max_up_dur7{d}(w) = nan;
%         end

        curDowns_w(curDowns_w==length(down_trans{d})) = [];
        curDowns_8(curDowns_8==length(down_trans8{d})) = [];
%         curDowns_5(curDowns_5==length(down_trans5{d})) = [];
%         curDowns_7(curDowns_7==length(down_trans7{d})) = [];

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
%         if ~isempty(curDowns_5)
%             max_down_dur5{d}(w) = max(down_state_dur5{d}(curDowns_5));
%         else
%             max_down_dur5{d}(w) = nan;
%         end
%         if ~isempty(curDowns_7)
%             max_down_dur7{d}(w) = max(down_state_dur7{d}(curDowns_7));
%         else
%             max_down_dur7{d}(w) = nan;
%         end


        avg_down_dur{d}(w) = mean(down_state_dur{d}(curDowns_w));
        avg_down_dur8{d}(w) = mean(down_state_dur8{d}(curDowns_8));
%         avg_down_dur5{d}(w) = mean(down_state_dur5{d}(curDowns_5));
%         avg_down_dur7{d}(w) = mean(down_state_dur7{d}(curDowns_7));

    end
    %

    des_times = find(up_state_dur8{d} > 5);
%     des_times = find(up_state_dur5{d} > 5);
    
    end_time = length(max_up_dur{d});

    plot(max_up_dur{d},'linewidth',2)
    hold on
    plot(max_up_dur8{d},'r','linewidth',2)
%     plot(max_up_dur5{d},'k','linewidth',2)
    xlim([0 end_time])
    plot(up_trans8{d}(des_times)/Fsd,ones(size(des_times))*5,'k*')
    cell_name = ['L' sess_data(d).layer '_' sess_data(d).cell_type '_' sess_data(d).region '_' sess_data(d).name];
    tname = ['C:\WC_Germany\Cortical_analysis\UDS_dur_raw\max_up_dur_new_method_005hz_' cell_name];
    print('-dpng',tname)
    close all


    %
    clear lf8 wcv_minus_spike 
    
end

save C:\WC_Germany\Cortical_analysis\UDS_dur_raw\UDS_raw_data up_state_dur* down_state_dur* up_trans* down_trans* 

