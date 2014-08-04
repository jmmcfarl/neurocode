clear all

load C:\WC_Germany\Persistent_activity\dir_tree_update
load C:\WC_Germany\Persistent_activity\UDS_dur_raw\UDS_raw_data
load C:\WC_Germany\Persistent_activity\sigmoid_fit\sig_fit_data
load C:\WC_Germany\Persistent_activity\lf8_period_f_data
load C:\WC_Germany\Persistent_activity\UDS_dist\synch_UDS_inds

niqf = 2016/2;
lcf = 0.05/niqf;
hcf = 2/niqf;
[b,a] = butter(2,[lcf hcf]);

dsf = 8;
Fsd = 2016/dsf;
backlag = 3*Fsd;
forwardlag = 10*Fsd
lags = -backlag:forwardlag;
% dlags = -backlag:backlag;

for d = 1:length(dir_array)
    cd(dir_array{d})
    pwd

    load used_data lf8 wcv_minus_spike

    lf8_f = filtfilt(b,a,lf8);
    wcv_f = filtfilt(b,a,wcv_minus_spike);

    lf8_d = downsample(lf8_f,dsf);
    wcv_d = downsample(wcv_f,dsf);

    lf8_d = zscore(lf8_d);
    wcv_d = zscore(wcv_d);

    wcv_up_ctrig_mat = zeros(length(synch_ups{d}),length(lags));
    lf8_up_ctrig_mat = wcv_up_ctrig_mat;

    new_up_mp = round(rlid_wup{d}/dsf);
    new_down_mp = round(rlid_wup{d}/dsf);
    new_up_lfp = round(rlid_lup{d}/dsf);
    
    
    %initialize
    secondary_lfp_states = cell(length(synch_ups{d}),1);
    mp_down_state = zeros(length(synch_ups{d}),1);
    lfp_pup = zeros(1,length(synch_ups{d}));
    lfp_pdown = lfp_pup;
    lfp_nup = lfp_pup;
    lfp_pupdur = lfp_pup;
    lfp_period_dur{d} = lfp_pup;

    for i = 1:length(synch_ups{d})

%         %find duration of MP up state in units of LFP periods
%         lfp_period_dur{d}(i) = lf8_period_f{d}(down_trans{d}(synch_ups{d}(i)))-lf8_period_f{d}(up_trans{d}(synch_ups{d}(i)));

%         %find first lfp down state after MP up transition
%         first_in_lfp_down = find(down_trans8{d} > up_trans{d}(synch_ups{d}(i)),1,'first');

%         if ~isempty(first_in_lfp_down)
%             %find all lfp up transitions that occur before mp comes back down
%             secondary_lfp_states{i} = find(down_trans8{d} < down_trans{d}(synch_ups{d}(i)) ...
%                 & up_trans8{d} > down_trans8{d}(first_in_lfp_down));
% 
%             %if there are secondary lfp up states, set up portions to nans
%             %for visualization purposes
%             if ~isempty(secondary_lfp_states{i})
%                 for q = 1:length(secondary_lfp_states{i})
%                     lf8_d(up_trans8{d}(secondary_lfp_states{i}(q)):down_trans8{d}(secondary_lfp_states{i}(q))) = nan;
%                 end
%             end
% 
%         end

%         %acquire MP up related lfp transitions
%         prev_lfp_up = find(up_trans8{d}<up_trans{d}(synch_ups{d}(i)),1,'last');
%         prev_lfp_down = find(down_trans8{d}<down_trans{d}(synch_ups{d}(i)),1,'last');
%         next_lfp_up = find(up_trans8{d}>down_trans{d}(synch_ups{d}(i)),1,'first');
%         if ~isempty(prev_lfp_up)
%             lfp_pup(i) = up_trans8{d}(prev_lfp_up);
%             lfp_pupdur(i) = up_state_dur8{d}(prev_lfp_up);
%         else
%             lfp_pup(i) = nan;
%             lfp_pupdur(i) = nan;
%         end
%         if ~isempty(prev_lfp_down)
%             lfp_pdown(i) = down_trans8{d}(prev_lfp_down);
%         else
%             lfp_pdown(i) = nan;
%         end
%         if ~isempty(next_lfp_up)
%             lfp_nup(i) = up_trans8{d}(next_lfp_up);
%         else
%             lfp_nup(i) = nan;
%         end
% 
%         %check if MP down transition occured during and LFP down state
%         if ~isempty(prev_lfp_down) & prev_lfp_down < length(up_trans8{d})
%             mp_down_state(i) = up_trans8{d}(prev_lfp_down+1) > down_trans{d}(synch_ups{d}(i));
%         else
%             mp_down_state(i) = nan;
%         end

        %calculate mup triggered averages
        if new_up_mp(synch_ups{d}(i)) > backlag & new_up_mp(synch_ups{d}(i)) ...
                < length(wcv_d)-forwardlag
            wcv_up_ctrig_mat(i,:) = wcv_d(new_up_mp(synch_ups{d}(i))-backlag: ...
                new_up_mp(synch_ups{d}(i))+forwardlag);
            lf8_up_ctrig_mat(i,:) = lf8_d(new_up_mp(synch_ups{d}(i))-backlag: ...
                new_up_mp(synch_ups{d}(i))+forwardlag);
        else
            wcv_up_ctrig_mat(i,:) = nan;
            lf8_up_ctrig_mat(i,:) = nan;
        end

    end

%     tsld = (down_trans{d}(synch_ups{d})-lfp_pdown)/Fsd; %time from mp down since last lfp down state
%     ttnu = (lfp_nup-down_trans{d}(synch_ups{d}))/Fsd; %time from mp down to next lfp up state
%     tslu = (up_trans{d}(synch_ups{d})-lfp_pup)/Fsd; %time from mp up since last lfp up
% 

    [dummy,up_order] = sort(up_state_dur{d}(synch_ups{d}));

    Fig = figure(1)
    clf
    set(Fig,'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
    set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
    %
    subplot(2,1,1)
    pcolor(lags/Fsd,(1:length(synch_ups{d}))/length(synch_ups{d}) ...
        ,wcv_up_ctrig_mat(up_order,:));shading flat;
    colorbar
    hold on
    plot(up_state_dur{d}(synch_ups{d}(up_order)),(1:length(synch_ups{d}))...
        /length(synch_ups{d}),'w','linewidth',2)
    line([0 0],[0 1],'Color','k')
    subplot(2,1,2)
    pcolor(lags/Fsd,(1:length(synch_ups{d}))/length(synch_ups{d})...
        ,lf8_up_ctrig_mat(up_order,:));shading flat;
    colorbar
    hold on
    plot(up_state_dur{d}(synch_ups{d}(up_order)),(1:length(synch_ups{d}))/length(synch_ups{d}),'w','linewidth',2)
    line([0 0],[0 1],'Color','k')
    tname = ['C:\WC_Germany\Persistent_activity\persistent_analysis\mp_up_trig_avg_sigfit' f_names{d}];
    print('-dpng',tname);
    close

%     %calculate fraction of MP up states with secondary states
%     sec_lfp_check{d} = zeros(length(synch_ups{d}),1);
%     num_sec_states{d} = zeros(length(synch_ups{d}),1);
% 
%     for i = 1:length(synch_ups{d})
% 
%         sec_lfp_check{d}(i) = ~isempty(secondary_lfp_states{i});
%         num_sec_states{d}(i) = length(secondary_lfp_states{i});
% 
%     end
% 
%     pers_fract(d) = nansum(sec_lfp_check{d})/length(sec_lfp_check{d});
%     mp_down_frac(d) = nansum(mp_down_state)/length(mp_down_state);
% 
% 
%     up_locking_cv(d) = nanstd(tslu)/nanmean(tslu);
%     down_locking_cv(d) = nanstd(tsld)/nanmean(tsld);
% 
%     hist(lfp_period_dur{d},100)
%     tname = ['C:\WC_Germany\Persistent_activity\persistent_analysis\lfp_fract_period_' f_names{d}];
%     print('-dpng',tname);
%     close

    clear near_lfp_ind
    clear prev_lfp_down

end

