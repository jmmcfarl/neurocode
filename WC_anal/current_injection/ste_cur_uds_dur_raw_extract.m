clear all

load C:\WC_Germany\current_injection\stellate\ste_current_dir
load C:\WC_Germany\current_injection\stellate\run_hist_thresh\UDS_dur_data

segDur = 20;
winSlide = 1;
dsf = 8;
Fsd = 2016/dsf;

thresh_smooth = 10;

niqf = 1/2;
[b,a] = butter(2,.005/niqf,'low');

for d = 1:length(dir_array)

    cd(dir_array{d})
    pwd

    load used_data lf8 wcv_minus_spike wcv
    
    sm_threshold_w{d} = filtfilt(b,a,threshold_w{d});
    sm_threshold_8{d} = filtfilt(b,a,threshold_8{d});
    
    disp('WCV')
    [up_state_dur{d},down_state_dur{d},up_trans{d},down_trans{d},sig_w{d},int_thresh_w{d},tot_chopped_up(d),tot_chopped_down(d)] = UDS_extract_run_hist(wcv_minus_spike,sm_threshold_w{d});
    disp('LFP')
    [up_state_dur8{d},down_state_dur8{d},up_trans8{d},down_trans8{d},sig_8{d},int_thresh_8{d},tot_chopped_up8(d),tot_chopped_down8(d)] = UDS_extract_run_hist(lf8,sm_threshold_8{d});

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

%     plot(max_up_dur{d},'linewidth',2)
%     hold on
%     plot(max_up_dur8{d},'r','linewidth',2)
%     xlim([0 end_time])
%     plot(up_trans8{d}(des_times)/Fsd,ones(size(des_times))*5,'k*')
%     tname = ['C:\WC_Germany\persistent_revised\UDS_dur_raw\max_up_dur_new_method_005hz_' f_names{d}];
%     print('-dpng',tname)
%     close all


    %
    clear lf8 wcv_minus_spike 
    
end

save C:\WC_Germany\current_injection\stellate\UDS_dur_raw\UDS_raw_data up_state_dur* down_state_dur* up_trans* down_trans* 

