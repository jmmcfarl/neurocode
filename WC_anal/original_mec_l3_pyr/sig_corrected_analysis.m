%% analyze precise timing relationships between LFP and MP state
%% transitions for each cell

clear all
load C:\WC_Germany\JMM_analysis_pyr\WCV_LFP_up_trans_sig_fit_data
load C:\WC_Germany\JMM_Analysis_pyr\dir_tree_update
load C:\WC_Germany\JMM_analysis_pyr\UDS_dur_run_hist_v2\up_per_data_10_28


Fs = 2016;
dsf = 4;
Fsd = Fs/dsf;

niqf = 2016/2;
lcf = 0.05/niqf;
hcf = 100/niqf;
[b,a] = butter(2,[lcf hcf]);

lags = -2*Fsd:2*Fsd;
maxlag = 2*Fsd;

lag_range = linspace(0,1,50);

for d = 1:length(dir_array)

    cd(dir_array{d});
    pwd

    load used_data wcv_minus_spike lf8 wcv

    wcv = filtfilt(b,a,wcv_minus_spike);
    wcv = downsample(wcv,dsf);
    wcv = zscore(wcv);
    wcv = wcv';

    lf8 = filtfilt(b,a,lf8);
    lf8 = downsample(lf8,dsf);
    lf8 = zscore(lf8);
    lf8 = lf8';

    up_trans{d} = up_trans{d}*2;
    up_trans8{d} = up_trans8{d}*2;
    down_trans{d} = down_trans{d}*2;
    down_trans8{d} = down_trans8{d}*2;

    t_90_wup{d} = round(t_90_wup{d}/4);
    t_10_wup{d} = round(t_10_wup{d}/4);

    t_90_lup{d} = round(t_90_lup{d}/4);
    t_10_lup{d} = round(t_10_lup{d}/4);

    %% find lfp up transitions when mp in down state
    lup_wdown_trans = [];
    lup_wdown_wtrans = [];
    for i = 1:length(up_trans8{d})

        %find preceding mp down transition
        prec_down = find(down_trans{d} < up_trans8{d}(i),1,'last');

        if ~isempty(prec_down) & prec_down ~= length(up_trans{d})

            if up_trans{d}(prec_down+1) > up_trans8{d}(i) & ~isnan(t_10_lup{d}(i))
                lup_wdown_trans = [lup_wdown_trans i];
                lup_wdown_wtrans = [lup_wdown_wtrans prec_down+1];
            end

        end


    end

    %% find lfp up transitions when mp is in up state
    lup_wup_trans = [];
    for i = 1:length(up_trans8{d})

        prec_up = find(up_trans{d} < up_trans8{d}(i),1,'last');
        if ~isempty(prec_up)
            if down_trans{d}(prec_up) > up_trans8{d}(i) & ~isnan(t_10_lup{d}(i))
                lup_wup_trans = [lup_wup_trans i];
            end
        end

    end

    %% calculate average 10% and 90% lags
    lag_10{d} = (t_10_wup{d}(lup_wdown_wtrans)-t_10_lup{d}(lup_wdown_trans))/Fsd;
    lag_90{d} = (t_90_wup{d}(lup_wdown_wtrans)-t_90_lup{d}(lup_wdown_trans))/Fsd;

    lag_10_dist(d,:) = hist(lag_10{d},lag_range);
    lag_90_dist(d,:) = hist(lag_90{d},lag_range);

    %% calculate LFP triggered MP averages
    wcv_lup_wdown_10 = zeros(length(lup_wdown_trans),length(lags));
    lfp_lup_wdown_10 = zeros(length(lup_wdown_trans),length(lags));
    wcv_lup_wdown_90 = zeros(length(lup_wdown_trans),length(lags));
    lfp_lup_wdown_90 = zeros(length(lup_wdown_trans),length(lags));

    for i = 1:length(lup_wdown_trans)
        if up_trans8{d}(lup_wdown_trans(i)) > maxlag & length(lf8)-up_trans8{d}(lup_wdown_trans(i)) > maxlag
            wcv_lup_wdown_10(i,:) = wcv(t_10_lup{d}(lup_wdown_trans(i))-maxlag:t_10_lup{d}(lup_wdown_trans(i))+maxlag);
            lfp_lup_wdown_10(i,:) = lf8(t_10_lup{d}(lup_wdown_trans(i))-maxlag:t_10_lup{d}(lup_wdown_trans(i))+maxlag);
            wcv_lup_wdown_90(i,:) = wcv(t_90_lup{d}(lup_wdown_trans(i))-maxlag:t_90_lup{d}(lup_wdown_trans(i))+maxlag);
            lfp_lup_wdown_90(i,:) = lf8(t_90_lup{d}(lup_wdown_trans(i))-maxlag:t_90_lup{d}(lup_wdown_trans(i))+maxlag);
        end
    end

    wcv_lup_wup_10 = zeros(length(lup_wup_trans),length(lags));
    lfp_lup_wup_10 = zeros(length(lup_wup_trans),length(lags));
    wcv_lup_wup_90 = zeros(length(lup_wup_trans),length(lags));
    lfp_lup_wup_90 = zeros(length(lup_wup_trans),length(lags));

    for i = 1:length(lup_wup_trans)
        if up_trans8{d}(lup_wup_trans(i)) > maxlag & length(lf8) - up_trans8{d}(lup_wup_trans(i)) > maxlag
            wcv_lup_wup_10(i,:) = wcv(t_10_lup{d}(lup_wup_trans(i))-maxlag:t_10_lup{d}(lup_wup_trans(i))+maxlag);
            lfp_lup_wup_10(i,:) = lf8(t_10_lup{d}(lup_wup_trans(i))-maxlag:t_10_lup{d}(lup_wup_trans(i))+maxlag);
            wcv_lup_wup_90(i,:) = wcv(t_90_lup{d}(lup_wup_trans(i))-maxlag:t_90_lup{d}(lup_wup_trans(i))+maxlag);
            lfp_lup_wup_90(i,:) = lf8(t_90_lup{d}(lup_wup_trans(i))-maxlag:t_90_lup{d}(lup_wup_trans(i))+maxlag);

        end
    end

%     subplot(2,1,1)
%     pcolor(lags,1:size(wcv_lup_wdown_10,1),wcv_lup_wdown_10);shading flat
%     subplot(2,1,2)
%     pcolor(lags,1:size(lfp_lup_wdown_10,1),lfp_lup_wdown_10);shading flat
%     t_names = ['C:\WC_Germany\JMM_analysis_pyr\sig_corrected_data\_lup_wdown_10' f_names{d}];
%     print('-dpng',t_names)
%     close all
% 
%     subplot(2,1,1)
%     pcolor(lags,1:size(wcv_lup_wdown_90,1),wcv_lup_wdown_90);shading flat
%     subplot(2,1,2)
%     pcolor(lags,1:size(lfp_lup_wdown_90,1),lfp_lup_wdown_90);shading flat
%     t_names = ['C:\WC_Germany\JMM_analysis_pyr\sig_corrected_data\_lup_wdown_90' f_names{d}];
%     print('-dpng',t_names)
%     close all
% 
%     subplot(2,1,1)
%     pcolor(lags,1:size(wcv_lup_wup_10,1),wcv_lup_wup_10);shading flat
%     subplot(2,1,2)
%     pcolor(lags,1:size(lfp_lup_wup_10,1),lfp_lup_wup_10);shading flat
%     t_names = ['C:\WC_Germany\JMM_analysis_pyr\sig_corrected_data\_lup_wup_10' f_names{d}];
%     print('-dpng',t_names)
%     close all
% 
%     
%         subplot(2,1,1)
%     pcolor(lags,1:size(wcv_lup_wup_90,1),wcv_lup_wup_90);shading flat
%     subplot(2,1,2)
%     pcolor(lags,1:size(lfp_lup_wup_90,1),lfp_lup_wup_90);shading flat
%     t_names = ['C:\WC_Germany\JMM_analysis_pyr\sig_corrected_data\_lup_wup_90' f_names{d}];
%     print('-dpng',t_names)
%     close all

    
end



