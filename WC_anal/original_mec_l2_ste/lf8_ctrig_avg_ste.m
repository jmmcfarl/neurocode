clear all

load C:\WC_Germany\JMM_analysis_ste\dir_tree_ste
load C:\WC_Germany\JMM_analysis_ste\UDS_dur_run_hist\data
load C:\WC_Germany\JMM_analysis_ste\lf8_period_f_data

niqf = 2016/2;
lcf = 0.05/niqf;
hcf = 2/niqf;
[b,a] = butter(2,[lcf hcf]);

dsf = 8;
Fsd = 2016/dsf;
backlag = 3*Fsd;
forwardlag = 10*Fsd
lags = -backlag:forwardlag;
dlags = -backlag:backlag;
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

    wcv_up_ctrig_mat = zeros(length(up_trans{d}),length(lags));
    lf8_up_ctrig_mat = wcv_up_ctrig_mat;

    %     near_lfp_ind = zeros(length(up_trans{d}),1);
    secondary_lfp_states = cell(length(up_trans{d}),1);
    mp_down_state = zeros(length(up_trans{d}),1);
lfp_pup = zeros(1,length(up_trans{d}));
lfp_pdown = lfp_pup;
lfp_nup = lfp_pup;
lfp_pupdur = lfp_pup;
lfp_period_dur{d} = lfp_pup;

    for i = 1:length(up_trans{d})

        %find closest lfp up transition
        %        [dummy,close_lf8up_ind] = min(abs(up_trans8{d}-up_trans{d}(i)));
        %        close_lf8up = up_trans8{d}(close_lf8up_ind);
        lfp_period_dur{d}(i) = lf8_period_f{d}(down_trans{d}(i))-lf8_period_f{d}(up_trans{d}(i));
        
        first_in_lfp_down = find(down_trans8{d} > up_trans{d}(i),1,'first');
        if ~isempty(first_in_lfp_down)
            %find all lfp up transitions that occur before mp comes back down
            secondary_lfp_states{i} = find(down_trans8{d} < down_trans{d}(i) & up_trans8{d} > down_trans8{d}(first_in_lfp_down));
%             if ~isempty(secondary_lfp_states{i})
%                 for q = 1:length(secondary_lfp_states{i})
%                     lf8_d(up_trans8{d}(secondary_lfp_states{i}(q)):down_trans8{d}(secondary_lfp_states{i}(q))) = nan;
%                 end
%             end
        end
        %does mp down occur during an lfp down
        prev_lfp_up = find(up_trans8{d}<up_trans{d}(i),1,'last');
        prev_lfp_down = find(down_trans8{d}<down_trans{d}(i),1,'last');
        next_lfp_up = find(up_trans8{d}>down_trans{d}(i),1,'first');
        if ~isempty(prev_lfp_up)
            lfp_pup(i) = up_trans8{d}(prev_lfp_up);
            lfp_pupdur(i) = up_state_dur8{d}(prev_lfp_up);
        else
            lfp_pup(i) = nan;
            lfp_pupdur(i) = nan;
        end
        if ~isempty(prev_lfp_down)
            lfp_pdown(i) = down_trans8{d}(prev_lfp_down);
        else
            lfp_pdown(i) = nan;
        end
        if ~isempty(next_lfp_up)
            lfp_nup(i) = up_trans8{d}(next_lfp_up);
        else
            lfp_nup(i) = nan;
        end
        
        if ~isempty(prev_lfp_down) & prev_lfp_down < length(up_trans8{d})
            mp_down_state(i) = up_trans8{d}(prev_lfp_down+1) > down_trans{d}(i);
        else
            mp_down_state(i) = nan;
        end
        if up_trans{d}(i) > backlag & up_trans{d}(i) < length(wcv_d)-forwardlag
            wcv_up_ctrig_mat(i,:) = wcv_d(up_trans{d}(i)-backlag:up_trans{d}(i)+forwardlag);
            lf8_up_ctrig_mat(i,:) = lf8_d(up_trans{d}(i)-backlag:up_trans{d}(i)+forwardlag);
        else
            wcv_up_ctrig_mat(i,:) = nan;
            lf8_up_ctrig_mat(i,:) = nan;
        end

    end

    %     wcv_down_ctrig_mat = zeros(length(down_trans{d}),length(dlags));
    %     lf8_down_ctrig_mat = wcv_down_ctrig_mat;
    %
    %         for i = 1:length(down_trans{d})
    %
    %
    % %             prev_lfp_down(i) = down_trans8{d}(find(down_trans8{d} < down_trans{d}(i),1,'first'));
    %
    %            if down_trans{d}(i) > forwardlag & down_trans{d}(i) < length(wcv_d)-backlag
    %                wcv_down_ctrig_mat(i,:) = wcv_d(down_trans{d}(i)-backlag:down_trans{d}(i)+backlag);
    %                lf8_down_ctrig_mat(i,:) = lf8_d(down_trans{d}(i)-backlag:down_trans{d}(i)+backlag);
    %            else
    %                wcv_down_ctrig_mat(i,:) = nan;
    %                lf8_down_ctrig_mat(i,:) = nan;
    %            end
    %
    %         end

    %     std_down_diff(d) = std(lfp_down_diff);
    %     std_up_diff(d) = std(lfp_up_diff);
    %    cv_down_diff(d) = std(lfp_down_diff)/mean(lfp_down_diff);
    %    cv_up_diff(d) = std(lfp_up_diff)/mean(lfp_up_diff);

    lfp_dur_calc = (lfp_pdown-lfp_pup)/Fsd;
    tsld = (down_trans{d}-lfp_pdown)/Fsd;
    ttnu = (lfp_nup-down_trans{d})/Fsd;
    tslu = (up_trans{d}-lfp_pup)/Fsd;
    
    up_dur_ratio{d} = up_state_dur{d}./lfp_dur_calc;
  

    [dummy,up_order] = sort(up_state_dur{d});
%     Fig = figure(1)
%     clf
%     set(Fig,'PaperUnits','centimeters');
%     set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
%     set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
%     %
%     subplot(2,1,1)
%     pcolor(lags/Fsd,(1:length(up_trans{d}))/length(up_trans{d}) ...
%         ,wcv_up_ctrig_mat(up_order,:));shading flat;
%     colorbar
%     hold on
%     plot(up_state_dur{d}(up_order),(1:length(up_trans{d}))/length(up_trans{d}),'w','linewidth',2)
%     line([0 0],[0 1],'Color','k')
%     subplot(2,1,2)
%     pcolor(lags/Fsd,(1:length(up_trans{d}))/length(up_trans{d})...
%         ,lf8_up_ctrig_mat(up_order,:));shading flat;
%     colorbar
%     hold on
%     plot(up_state_dur{d}(up_order),(1:length(up_trans{d}))/length(up_trans{d}),'w','linewidth',2)
%     line([0 0],[0 1],'Color','k')
%     %     for i = 1:length(up_trans{d})
%     %        if ~isempty(secondary_lfp_states{i})
%     %           plot((up_trans8{d}(secondary_lfp_states{up_order(i)})-up_trans{d}(up_order(i)))/Fsd,ones(size(secondary_lfp_states{up_order(i)}))*i,'k.')
%     %           plot((down_trans8{d}(secondary_lfp_states{up_order(i)})-up_trans{d}(up_order(i)))/Fsd,ones(size(secondary_lfp_states{up_order(i)}))*i,'w.')
%     %        end
%     %     end
%     tname = ['C:\WC_Germany\JMM_Analysis_ste\lf8_up_ctrig\' f_names{d}];
%     print('-dpng',tname);
%     close

    %calculate fraction of MP up states with secondary states
    sec_lfp_check{d} = zeros(length(up_trans{d}),1);
    num_sec_states{d} = zeros(length(up_trans{d}),1);

    for i = 1:length(up_trans{d})

        sec_lfp_check{d}(i) = ~isempty(secondary_lfp_states{i});
        num_sec_states{d}(i) = length(secondary_lfp_states{i});

    end

    pers_fract(d) = nansum(sec_lfp_check{d})/length(sec_lfp_check{d});
    mp_down_frac(d) = nansum(mp_down_state)/length(mp_down_state);

hrange = linspace(0,3,50);

% no_pers_states = find(num_sec_states{d}==0);
% pers_states = find(num_sec_states{d} > 0);
% plot(up_state_dur{d}(no_pers_states),lfp_pupdur(no_pers_states),'o')
% hold on
% plot(up_state_dur{d}(pers_states),lfp_pupdur(pers_states),'ro')
% ylim([0 3])
% tname = ['C:\WC_Germany\JMM_Analysis_Pyr\lf8_up_ctrig\up_dur_cor_PS_' f_names{d}];
% print('-dpng',tname);
% close
% 
% plot(tslu(no_pers_states),tsld(no_pers_states),'o')
% hold on
% plot(tslu(pers_states),tsld(pers_states),'ro')
% xlim([0 3])
% ylim([0 3])
% tname = ['C:\WC_Germany\JMM_Analysis_Pyr\lf8_up_ctrig\tsld_tslu_PS_'
% f_names{d}];
% print('-dpng',tname);
% close

up_locking_cv(d) = nanstd(tslu)/nanmean(tslu);
down_locking_cv(d) = nanstd(tsld)/nanmean(tsld);

lfp_period_dur{d}(lfp_period_dur{d}<0) = [];
hist(lfp_period_dur{d},100)
tname = ['C:\WC_Germany\JMM_Analysis_ste\lf8_up_ctrig\lfp_fract_period_' f_names{d}];
print('-dpng',tname);
close


    d
    %
    %   figure
    %       subplot(2,1,1)
    %     pcolor(dlags/Fsd,1:length(down_trans{d}),wcv_down_ctrig_mat(up_order,:));shading flat
    %     subplot(2,1,2)
    %     pcolor(dlags/Fsd,1:length(down_trans{d}),lf8_down_ctrig_mat(up_order,:));shading flat
    %     tname = ['C:\WC_Germany\JMM_Analysis_Pyr\lf8_up_ctrig\down_' f_names{d}];
    %     print('-dpng',tname);
    %     close
    %

    clear near_lfp_ind
    clear prev_lfp_down

end

