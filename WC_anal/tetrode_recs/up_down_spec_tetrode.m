clear all
cd C:\WC_Germany\Tetrode_recs\2008-12-16
load UDS_dur avg* up_state_dur* down_state_dur* up_trans* down_trans* 

dsf = 2;
TW = 1;
params.Fs = 2016/dsf;
params.err = [1 0.01];
params.fpass = [20 200];
params.tapers = [TW 2*TW-1];
bad_f = [50 100 150 200];
T = 0.5;
W = TW/T;
% W = 0.04;
win = 100;
niqf = 2016/2;
lcf = 10/niqf;
hcf = 250/niqf;
[b,a] = butter(2,[lcf hcf]);

movingwin = [T T];
d = 1

    load lfp_data
    
    wcv_minus_spike = lf1_data;
    lf8 = lf12_data;

    wcv_f = filtfilt(b,a,wcv_minus_spike);
    lf8_f = filtfilt(b,a,lf8);
    down_w = downsample(wcv_f,dsf);
    down_8 = downsample(lf8_f,dsf);
%     down_w = downsample(wcv_minus_spike,dsf);
%     down_8 = downsample(lf8,dsf);

    down_w = zscore(down_w);
    down_8 = zscore(down_8);

    up_trans8 = up_trans8*4;
    down_trans8 = down_trans8*4;
    up_trans = up_trans*4;
    down_trans = down_trans*4;

    up_markers8 = [up_trans8' down_trans8'];
    up_markersw = [up_trans' down_trans'];
    down_markers8 = [down_trans8(1:end-1)' up_trans8(2:end)'];
    down_markersw = [down_trans(1:end-1)' down_trans(2:end)'];
    
    down_w = down_w';
down_8 = down_8';

    [S8_lup,f]= mtspectrumc_unequal_length_trials(down_8, movingwin, params, up_markers8);
    [S8_ldown,f]= mtspectrumc_unequal_length_trials(down_8, movingwin, params, down_markers8);
        [S8_up,f]= mtspectrumc_unequal_length_trials(down_8, movingwin, params, up_markersw);
    [S8_down,f]= mtspectrumc_unequal_length_trials(down_8, movingwin, params, down_markersw);
    [Sw_up,f]= mtspectrumc_unequal_length_trials(down_w, movingwin, params, up_markersw);
    [Sw_down,f]= mtspectrumc_unequal_length_trials(down_w, movingwin, params, down_markersw);
        [Sw_lup,f]= mtspectrumc_unequal_length_trials(down_w, movingwin, params, up_markers8);
    [Sw_ldown,f]= mtspectrumc_unequal_length_trials(down_w, movingwin, params, down_markers8);

%     [Cmn_up{d},Phimn_up{d},Smn,Smm,f,ConfC_up{d},PhiStd_up{d},Cerr_up{d}] = coherencyc_unequal_length_trials([down_8 down_w], movingwin, params,up_markersw );
%     [Cmn_down{d},Phimn_down{d},Smn,Smm,f,ConfC_down{d},PhiStd_down{d},Cerr_down{d}] = coherencyc_unequal_length_trials([down_8 down_w], movingwin, params,down_markersw );
%     [Cmn_lup{d},Phimn_lup{d},Smn,Smm,f,ConfC_lup{d},PhiStd_lup{d},Cerr_lup{d}] = coherencyc_unequal_length_trials([down_8 down_w], movingwin, params,up_markers8 );
%     [Cmn_ldown{d},Phimn_ldown{d},Smn,Smm,f,ConfC_ldown{d},PhiStd_ldown{d},Cerr_ldown{d}] = coherencyc_unequal_length_trials([down_8 down_w], movingwin, params,down_markers8 );

    exclude_f = [];
    for i = 1:length(bad_f)
        exclude_f = [exclude_f find(abs(f-bad_f(i)) < 2*W)];
    end
    S8_lup(exclude_f) = nan;
    S8_ldown(exclude_f) = nan;
    S8_up(exclude_f) = nan;
    S8_down(exclude_f) = nan;
    Sw_up(exclude_f) = nan;
    Sw_down(exclude_f) = nan;
        Sw_lup(exclude_f) = nan;
    Sw_ldown(exclude_f) = nan;

% Cmn_up{d}(exclude_f) = nan;
% Cmn_down{d}(exclude_f) = nan;
% Cmn_lup{d}(exclude_f) = nan;
% Cmn_ldown{d}(exclude_f) = nan;

% %     %% plot and save results
% subplot(2,1,1)
% plot(f,Cmn_up{d},'linewidth',2)
% hold on
% plot(f,Cmn_down{d},'r','linewidth',2)
% legend('MP UP state','MP Down State')
% line([0 max(f)],[ConfC_up{d} ConfC_up{d}],'Color','b')
% line([0 max(f)],[ConfC_down{d} ConfC_down{d}],'Color','r')
% subplot(2,1,2)
% plot(f,Cmn_lup{d},'linewidth',2)
% hold on
% plot(f,Cmn_down{d},'r','linewidth',2)
% legend('LFP UP state','LFP Down state')
% line([0 max(f)],[ConfC_lup{d} ConfC_lup{d}],'Color','b')
% line([0 max(f)],[ConfC_ldown{d} ConfC_ldown{d}],'Color','r')
% t_names = ['C:\WC_Germany\JMM_analysis_pyr\up_down_spectra\coh_20_400_' f_names{d}];
% print('-dpng',t_names)
% close all

%     subplot(2,1,1)
%     plot(f,10*log10(S8_up{d}),'linewidth',2)
%     hold on
%     plot(f,10*log10(S8_down{d}),'r','linewidth',2)
%     legend('Up','Down')
%     %     plot(f,10*log10(Serr8_up{d}(1,:)),'--')
%     %     plot(f,10*log10(Serr8_up{d}(2,:)),'--')
%     %     plot(f,10*log10(Serr8_down{d}(1,:)),'r--')
%     %     plot(f,10*log10(Serr8_down{d}(2,:)),'r--')
%     title('LFP')
%     subplot(2,1,1)
%     plot(f,10*log10(Sw_up{d}),'linewidth',2)
%     hold on
%     plot(f,10*log10(Sw_down{d}),'r','linewidth',2)
%     legend('Up','Down')
%     subplot(2,1,2)
%     plot(f,10*log10(Sw_up{d})-10*log10(Sw_down{d}),'k','linewidth',2)
    
%     %     plot(f,10*log10(Serrw_up{d}(1,:)),'--')
%     %     plot(f,10*log10(Serrw_up{d}(2,:)),'--')
%     %     plot(f,10*log10(Serrw_down{d}(1,:)),'r--')
%     %     plot(f,10*log10(Serrw_down{d}(2,:)),'r--')
%     title('MP')
%     t_names = ['C:\WC_Germany\JMM_analysis_pyr\up_down_spectra\up_down_pow_5_400' f_names{d}];
%     print('-dpng',t_names)
%     close all

    S8_diff = 10*log10(S8_up)-10*log10(S8_down);
    Sw_diff = 10*log10(Sw_up)-10*log10(Sw_down);

%     subplot(2,1,1)
%     plot(f,10*log10(S8_up{d})-10*log10(S8_down{d}),'linewidth',2)
%     hold on
%     plot(f,10*log10(Sw_up{d})-10*log10(Sw_down{d}),'r','linewidth',2)
%     legend('LFP','MP')
%     title('Up-Down power')
%     subplot(2,1,2)
%     plot(f,S8_diff{d}-Sw_diff{d},'k','linewidth',2)
%     t_names = ['C:\WC_Germany\JMM_analysis_pyr\up_down_spectra\diff_20_100' f_names{d}];
%     print('-dpng',t_names)
%     close all

    mean_hf_diff(d) = nanmean(S8_diff-Sw_diff);
mean_lup8_pow(d) = nanmean(10*log10(S8_lup));
mean_ldown8_pow(d) = nanmean(10*log10(S8_ldown));
mean_ldiff8_pow(d) = nanmean(10*log10(S8_lup)-10*log10(S8_ldown));
mean_up8_pow(d) = nanmean(10*log10(S8_up));
mean_down8_pow(d) = nanmean(10*log10(S8_down));
mean_diff8_pow(d) = nanmean(10*log10(S8_up)-10*log10(S8_down));
    mean_up_pow(d) = nanmean(10*log10(Sw_up));
    mean_down_pow(d) = nanmean(10*log10(Sw_down));
    mean_diff_pow(d) = nanmean(10*log10(Sw_up)-10*log10(Sw_down));
    mean_lup_pow(d) = nanmean(10*log10(Sw_lup));
    mean_ldown_pow(d) = nanmean(10*log10(Sw_ldown));
    mean_ldiff_pow(d) = nanmean(10*log10(Sw_lup)-10*log10(Sw_ldown));
% mean_coh_up(d) = nanmean(Cmn_up{d});
% mean_coh_down(d) = nanmean(Cmn_down{d});
% mean_coh_diff(d) = nanmean(Cmn_up{d}-Cmn_down{d});

%% calculate average power in different bands

save up_down_spectra mean*

