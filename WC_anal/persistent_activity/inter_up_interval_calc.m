clear all
close all
load C:\WC_Germany\Persistent_activity\dir_tree_update
load C:\WC_Germany\Persistent_activity\UDS_dur_raw\UDS_raw_data
load C:\WC_Germany\Persistent_activity\UDS_dist\synch_UDS_inds
load('C:\WC_Germany\Persistent_activity\persistent_analysis\quant_data_12_05.mat')
load('C:\WC_Germany\Persistent_activity\persistent_analysis\quant_data_12_05.mat')
niqf = 2016/2;
dsf = 8;
Fsd = 2016/8;
lcf = 0.05/niqf;
hcf = 4/niqf;
[b,a] = butter(2,[lcf hcf]);
for d = 1:length(dir_array)
% cd(dir_array{d})
% pwd
%     load used_data lf8
%     
%     lf8_f = filtfilt(b,a,lf8);
%     lf8_d = downsample(lf8_f,dsf);
%     lf8_d = zscore(lf8_d);
    
    
    iui{d} = zeros(size(synch_ups{d}));
    near_l_dur{d} = zeros(size(synch_ups{d}));
    near_ld_dur{d} = zeros(size(synch_ups{d}));
    tslu{d} = zeros(size(synch_ups{d}));
    extra_pers{d}= zeros(size(synch_ups{d}));
    next_peak{d} = zeros(size(synch_ups{d}));
    ave_iui{d} = zeros(size(synch_ups{d}));
    
    for i = 1:length(synch_ups{d})

               prev_lfp_up = up_trans8{d}(find(up_trans8{d} < up_trans{d}(synch_ups{d}(i)),1,'last'));
%                if ~isempty(prev_lfp_up)
%                next_peak{d}(i) = findpeaks(lf8_d(prev_lfp_up:end),'npeaks',1);
%                else
%                    next_peak{d}(i) = nan;
%                end
        [dummy,near_lfp_up] = min(abs(up_trans8{d}-up_trans{d}(synch_ups{d}(i))));
        
        cur_lfp_ups = near_lfp_up-1+find(up_trans8{d}(near_lfp_up:end)<down_trans{d}(synch_ups{d}(i)));
        if length(cur_lfp_ups) > 1
            ave_iui{d}(i) = mean(diff(up_trans8{d}(cur_lfp_ups)))/Fsd;
        elseif length(cur_lfp_ups) ==1 & near_lfp_up < length(up_trans8{d})-1
            ave_iui{d}(i) = (up_trans8{d}(near_lfp_up+1)-up_trans8{d}(near_lfp_up))/Fsd;;
        else
            ave_iui{d}(i) = nan;
        end
        
        
                    near_lfp_down = find(down_trans8{d}<down_trans{d}(synch_ups{d}(i)),1,'last');
        if ~isempty(near_lfp_up) & near_lfp_up < length(up_trans8{d})-1
            iui{d}(i) = (up_trans8{d}(near_lfp_up+1)-up_trans8{d}(near_lfp_up))/Fsd;
            near_l_dur{d}(i) = up_state_dur8{d}(near_lfp_up);
            near_ld_dur{d}(i) = down_state_dur8{d}(near_lfp_up);
            tslu{d}(i) = (up_trans{d}(synch_ups{d}(i))-up_trans8{d}(near_lfp_up))/Fsd;
        else
            iui{d}(i) = nan;
            near_l_dur{d}(i) = nan;
            tslu{d}(i) = nan;
            near_ld_dur{d}(i) = nan;
        end
        if ~isempty(near_lfp_up) & near_lfp_up < length(up_trans8{d})-2
            iui2{d}(i) = (up_trans8{d}(near_lfp_up+2)-up_trans8{d}(near_lfp_up))/Fsd;
        else
            iui2{d}(i) = nan;
        end

        if ~isempty(near_lfp_up) & ~isempty(near_lfp_down)
            extra_pers{d}(i) = (down_trans{d}(synch_ups{d}(i))-down_trans8{d}(near_lfp_down))/Fsd;
        else
            extra_pers{d}(i) = nan;
        end
        
    end

    cur_lfp_dur = lfp_period_dur{d};
    nps = find(cur_lfp_dur < 1);
    pas = find(cur_lfp_dur > 1);
    ones = find(cur_lfp_dur > 1 & cur_lfp_dur < 2);
    twos = find(cur_lfp_dur > 2 & cur_lfp_dur < 3);
    vlong = find(cur_lfp_dur > 3);

    mean_iui_nps(d) = nanmean(iui{d}(nps));
    cv_iui_nps(d) = mean_iui_nps(d)/nanstd(iui{d}(nps));
    mean_iui_pas(d) = nanmean(iui{d}(pas));
    cv_iui_pas(d) = nanstd(iui{d}(pas))/mean_iui_pas(d);
    mean_iui_ones(d) = nanmean(iui{d}(ones));
    cv_iui_ones(d) = mean_iui_ones(d)/nanstd(iui{d}(ones));
    mean_iui_twos(d) = nanmean(iui{d}(twos));
    cv_iui_twos(d) = mean_iui_twos(d)/nanstd(iui{d}(twos));
    mean_iui_vlong(d) = nanmean(iui{d}(vlong));
    cv_iui_vlong(d) = mean_iui_vlong(d)/nanstd(iui{d}(vlong));

    cv_iui_nps(d) = 1/cv_iui_nps(d);
    cv_iui_ones(d) = 1/cv_iui_ones(d);
    cv_iui_twos(d) = 1/cv_iui_twos(d);
    cv_iui_vlong(d) = 1/cv_iui_vlong(d);

%     figure
%     plot(near_ld_dur{d}(nps)-near_l_dur{d}(nps),up_state_dur{d}(synch_ups{d}(nps)),'.')
%     hold on
%     plot(near_ld_dur{d}(ones)-near_l_dur{d}(ones),up_state_dur{d}(synch_ups{d}(ones)),'r.')
%     plot(near_ld_dur{d}(twos)-near_l_dur{d}(twos),up_state_dur{d}(synch_ups{d}(twos)),'g.')
%     plot(near_ld_dur{d}(vlong)-near_l_dur{d}(vlong),up_state_dur{d}(synch_ups{d}(vlong)),'m.')
% pause
% close


%     figure
%     plot(next_peak{d}(nps),tslu{d}(nps),'.')
%     hold on
%     plot(next_peak{d}(ones),tslu{d}(ones),'r.')
% %     plot(next_peak{d}(twos),up_state_dur{d}(synch_ups{d}(twos)),'g.')
% %     plot(next_peak{d}(vlong),up_state_dur{d}(synch_ups{d}(vlong)),'m.')
% pause
% close


%     figure
%     plot(iui{d}(nps),up_state_dur{d}(synch_ups{d}(nps)),'.')
%     hold on
%     plot(iui{d}(ones),up_state_dur{d}(synch_ups{d}(ones)),'r.')
%     plot(iui{d}(twos),up_state_dur{d}(synch_ups{d}(twos)),'g.')
%     plot(iui{d}(vlong),up_state_dur{d}(synch_ups{d}(vlong)),'m.')
%     legend('No persistence','One LFP Cycle','Two LFP Cycles','More than 2')
%     xlabel('Inter-up-interval (s)')
%     ylabel('MP up state duration(s)')
%     d
%     pause
%     close
          figure
        plot(cur_lfp_dur(nps),ave_iui{d}(nps),'.')
        hold on
        plot(cur_lfp_dur(ones),ave_iui{d}(ones),'r.')
        plot(cur_lfp_dur(twos),ave_iui{d}(twos),'g.')
        plot(cur_lfp_dur(vlong),ave_iui{d}(vlong),'m.')
      pause
      close

%           figure
%         plot(near_l_dur{d}(nps),extra_pers{d}(nps),'.')
%         hold on
%         plot(near_l_dur{d}(ones),extra_pers{d}(ones),'r.')
%         plot(near_l_dur{d}(twos),extra_pers{d}(twos),'g.')
%         plot(near_l_dur{d}(vlong),extra_pers{d}(vlong),'m.')
%       pause
%       close

%           figure
%         plot(extra_pers{d}(nps),up_state_dur{d}(synch_ups{d}(nps)),'.')
%         hold on
%         plot(extra_pers{d}(ones),up_state_dur{d}(synch_ups{d}(ones)),'r.')
%         plot(extra_pers{d}(twos),up_state_dur{d}(synch_ups{d}(twos)),'g.')
%         plot(extra_pers{d}(vlong),up_state_dur{d}(synch_ups{d}(vlong)),'m.')
%       pause
%       close

%             plot(near_l_dur{d}(nps),iui{d}(nps),'.')
%         hold on
%         plot(near_l_dur{d}(ones),iui{d}(ones),'r.')
%         plot(near_l_dur{d}(twos),iui{d}(twos),'g.')
%         plot(near_l_dur{d}(vlong),iui{d}(vlong),'m.')
% 
%         pause
%     close

%     scatter3(iui{d}(nps),up_state_dur{d}(synch_ups{d}(nps)),near_l_dur{d}(nps),'.')
%     hold on
%     scatter3(iui{d}(ones),up_state_dur{d}(synch_ups{d}(ones)),near_l_dur{d}(ones),'r.')
%     scatter3(iui{d}(twos),up_state_dur{d}(synch_ups{d}(twos)),near_l_dur{d}(twos),'g.')
%     scatter3(iui{d}(vlong),up_state_dur{d}(synch_ups{d}(vlong)),near_l_dur{d}(vlong),'m.')

end