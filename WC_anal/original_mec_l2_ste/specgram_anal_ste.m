clear all

load C:\WC_Germany\JMM_Analysis_ste\dir_tree_ste

dsf = 8;
params.Fs = 2016/dsf;
params.tapers = [4 7];
params.err = [2 .05];
params.fpass = [0 10];
window = [60 5];

niqf = 2016/2;
% lcf = .01/niqf;

% [b,a] = butter(2,lcf,'high');

for d = 1:length(dir_array)
    
    disp(sprintf('session %d',d))
    cd(dir_array{d});
    pwd
    
    load used_data wcv_minus_spike lf8
    
%     wcv_minus_spike = filtfilt(b,a,wcv_minus_spike);
%     lf8 = filtfilt(b,a,lf8);
    
    down_w = downsample(wcv_minus_spike,dsf);
    down_8 = downsample(lf8,dsf);
    down_w = zscore(down_w);
    down_8 = zscore(down_8);
    
    [Sw{d},t{d},f{d}] = mtspecgramc(down_w,window,params);
    [S8{d},t{d},f{d}] = mtspecgramc(down_8,window,params);
    [C{d},phi{d},S12,S1,S2,tp,fp,confC,phistd{d},Cerr] = cohgramc(down_w,down_8,window,params);
    
%     max_freq = find(f{d}>1,1,'first');
%     min_freq_w = find(f{d}>0.05,1,'first');
%     min_freq_8 = find(f{d}>0.1,1,'first');
%     
%     [max_val_C{d},max_freq_C{d}] = max(C{d}(:,min_freq_w:max_freq),[],2);
%     [max_val_Sw{d},max_freq_Sw{d}] = max(Sw{d}(:,min_freq_w:max_freq),[],2);
%     [max_val_S8{d},max_freq_S8{d}] = max(S8{d}(:,min_freq_8:max_freq),[],2);
%     max_freq_C{d} = max_freq_C{d} + min_freq_w;
%     max_freq_Sw{d} = max_freq_Sw{d}+min_freq_w;
%     max_freq_S8{d} = max_freq_S8{d}+min_freq_8;
   
    
    Fig = figure(1)
    clf
    set(Fig,'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [40 40]);% paper size is in [width height] format
    set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
    subplot(3,1,1)
    pcolor(t{d},f{d},10*log10(Sw{d}'));shading flat
%     hold on
%     plot(t{d},jmm_smooth_1d(f{d}(max_freq_Sw{d}),4),'k','linewidth',2)
    ylim([0 10])
    caxis([-30 6]); 
    title('WCV')
    subplot(3,1,2)
    pcolor(t{d},f{d},10*log10(S8{d}'));shading flat
%     hold on
%     plot(t{d},jmm_smooth_1d(f{d}(max_freq_S8{d}),4),'k','linewidth',2)
    ylim([0 10])
    caxis([-30 6]);
    title('LFP')
    subplot(3,1,3)
    pcolor(t{d},f{d},C{d}');shading flat
%     hold on
%     plot(t{d},jmm_smooth_1d(f{d}(max_freq_C{d}),4),'k','linewidth',2)
    ylim([0 10])
   
%     subplot(4,1,4)
%     plot(max_up_dur{d},'linewidth',2)
%     hold on
%     plot(max_up_dur8{d},'r','linewidth',2)
    tname = ['C:\WC_Germany\JMM_Analysis_ste\specgram\wide_band' f_names{d}];
    print('-dpng',tname);
    close

%plot coherency and phase
%     Fig = figure(1)
%     clf
%     set(Fig,'PaperUnits','centimeters');
%     set(gcf, 'PaperSize', [40 40]);% paper size is in [width height] format
%     set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
%     subplot(2,1,1)
%     pcolor(t{d},f{d},C{d}');shading flat
%     hold on
%     colorbar
%     plot(t{d},jmm_smooth_1d(f{d}(max_freq_C{d}),2),'k','linewidth',2)
%     ylim([0 2])
%     subplot(2,1,2)
%     pcolor(t{d},f{d},360/2/pi*phi{d}');shading flat
%     hold on
%     plot(t{d},jmm_smooth_1d(f{d}(max_freq_C{d}),2),'k','linewidth',2)
%     ylim([0 2]);
%     colorbar
%     tname = ['C:\WC_Germany\JMM_Analysis_pyr\spectra\phase' f_names{d}];
%     print('-dpng',tname);
%     close
% 
% 
% 
    save C:\WC_Germany\JMM_Analysis_ste\specgram\specgram_data_60 C S8 Sw t f phi
    
end