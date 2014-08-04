clear all
close all
cd C:\WC_Germany\current_injection\Simultaneous_LFP

part_rec{1} = 'C:\WC_Germany\april_09_data\2009-04-05_CWC_LFP\2009-4-5_18-26-31';
part_rec{2} = 'C:\WC_Germany\april_09_data\2009-04-07\2009-4-7-19';
part_rec{3} = 'C:\WC_Germany\april_09_data\2009-04-13_A\2009-4-13-18';
part_rec{4} = 'C:\WC_Germany\april_09_data\2009-04-21_CWC_LFP\2009-4-21_22-20-11';

part_name{1} = '2009-4-5';
part_name{2} = '2009-4-7';
part_name{3} = '2009-4-13';
part_name{4} = '2009-4-21';


Fs = 2016;
niqf = Fs/2;
dsf = 8;
Fsd = Fs/dsf;
lcf = 0.2/niqf;
hcf = 40/niqf;

[b,a] = butter(2,[lcf hcf]);

window = round(10*Fsd);
noverlap = window/2;

for d = 1:length(part_rec)
    
    cd(part_rec{d})
    pwd
    load raw_data
    wcv_d = filtfilt(b,a,wcv_minus_spike);
    wcv_d = zscore(downsample(wcv_d,dsf));
    lf8_d = filtfilt(b,a,lf8);
    lf8_d = zscore(downsample(lf8_d,dsf));
    lf3_d = filtfilt(b,a,lf3);
    lf3_d = zscore(downsample(lf3_d,dsf));
    
    t = (1:length(wcv_d))/Fsd;
    
    [Sw,F,T,Pw] = spectrogram(wcv_d,window,noverlap,[],Fsd);
    [S8,F,T,P8] = spectrogram(lf8_d,window,noverlap,[],Fsd);
    [S3,F,T,P3] = spectrogram(lf3_d,window,noverlap,[],Fsd);

    load Events_File
    load raw_sync_times
    synct_d = downsample(synct,dsf);
    plot(synct_d,wcv_d);hold on
    
        stim_on = find(Events_Nttls == -1);
plot(Events_TimeStamps(stim_on),ones(size(stim_on)),'k.')
    cd C:\WC_Germany\current_injection\Simultaneous_LFP

    
    
%     Fig=figure(1);
%     set(Fig,'PaperUnits','centimeters');
%     set(gcf, 'PaperSize', [40 20]);
%     set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
%     subplot(2,1,1), title('MP')
%     pcolor(T,F,10*log10(Pw));shading flat, ylim([0 10]), caxis([-60 5])
%     subplot(2,1,2), title('LF8')
%     pcolor(T,F,10*log10(P8));shading flat, ylim([0 10]), caxis([-60 5])
%     t_name = ['liny_specgram8_compare_' part_name{d}];
%     print('-dpng',t_name)
%     subplot(2,1,1), set(gca,'yscale','log'), ylim([0.2 10])
%     subplot(2,1,2), set(gca,'yscale','log'), ylim([0.2 10])
%     t_name = ['logy_specgram8_compare_' part_name{d}];
%     print('-dpng',t_name), close
%  
%     Fig=figure(1);
%     set(Fig,'PaperUnits','centimeters');
%     set(gcf, 'PaperSize', [40 20]);
%     set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
%     subplot(2,1,1), title('MP')
%     pcolor(T,F,10*log10(Pw));shading flat, ylim([0 10]), caxis([-60 5])
%     subplot(2,1,2), title('LF8')
%     pcolor(T,F,10*log10(P3));shading flat, ylim([0 10]), caxis([-60 5])
%     t_name = ['liny_specgram3_compare_' part_name{d}];
%     print('-dpng',t_name)
%     subplot(2,1,1), set(gca,'yscale','log'), ylim([0.2 10])
%     subplot(2,1,2), set(gca,'yscale','log'), ylim([0.2 10])
%     t_name = ['logy_specgram3_compare_' part_name{d}];
%     print('-dpng',t_name), close
    
end
