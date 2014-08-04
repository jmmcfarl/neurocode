clear all
close all
cd C:\WC_Germany\current_injection\stellate
load C:\WC_Germany\JMM_analysis_ste\dir_tree_ste
dir_array2 = dir_array;
dir_array2(6) = [];
dir_array2 = dir_array2(4:15);
load ste_current_dir

Fs = 2016;
niqf = Fs/2;
dsf = 8;
Fsd = Fs/dsf;
lcf = 0.2/niqf;
hcf = 40/niqf;

[b,a] = butter(2,[lcf hcf]);

window = round(10*Fsd);
noverlap = window/2;

for d = 1:length(dir_array)
    
        cd(dir_array2{d})
    load used_data lf8 lf3
    lf8_d = filtfilt(b,a,lf8);
    m8 = mean(lf8_d);
    s8 = std(lf8_d);
    lf8_d = zscore(downsample(lf8_d,dsf));
    lf3_d = filtfilt(b,a,lf3);
    m3 = mean(lf3_d);
    s3 = std(lf3_d);
    lf3_d = zscore(downsample(lf3_d,dsf));
    [S8,F2,T2,P82] = spectrogram(lf8_d,window,noverlap,[],Fsd);
    [S3,F2,T2,P32] = spectrogram(lf3_d,window,noverlap,[],Fsd);
    clear lf8* lf3*
    cd(dir_array{d})
    pwd
    load used_data lf8 lf3 wcv_minus_spike
    wcv_d = filtfilt(b,a,wcv_minus_spike);
    wcv_d = zscore(downsample(wcv_d,dsf));
    lf8_d = filtfilt(b,a,lf8);
    lf8_d = downsample(lf8_d,dsf);
    lf3_d = filtfilt(b,a,lf3);
    lf3_d = downsample(lf3_d,dsf);
    lf8_d = (lf8_d - m8)/s8;
    lf3_d = (lf3_d - m3)/s3;
    
    t = (1:length(wcv_d))/Fsd;
    
    [Sw,F,T,Pw] = spectrogram(wcv_d,window,noverlap,[],Fsd);
    [S8,F,T,P8] = spectrogram(lf8_d,window,noverlap,[],Fsd);
    [S3,F,T,P3] = spectrogram(lf3_d,window,noverlap,[],Fsd);


    
    cd C:\WC_Germany\current_injection\stellate

    
    
    Fig=figure(1);
    set(Fig,'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [40 20]);
    set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
    subplot(3,1,1), title('MP')
    pcolor(T,F,10*log10(Pw));shading flat, ylim([0 10]), caxis([-60 5])
    subplot(3,1,2), title('LF8')
    pcolor(T,F,10*log10(P8));shading flat, ylim([0 10]), caxis([-60 5])
    subplot(3,1,3), title('LF8 Prior')
    pcolor(T2,F2,10*log10(P82));shading flat, ylim([0 10]), caxis([-60 5])  
    t_name = ['liny_specgram8_compare_' f_names{d}];
    print('-dpng',t_name)
    subplot(3,1,1), set(gca,'yscale','log'), ylim([0.2 10])
    subplot(3,1,2), set(gca,'yscale','log'), ylim([0.2 10])
    subplot(3,1,3), set(gca,'yscale','log'), ylim([0.2 10])
    t_name = ['logy_specgram8_compare_' f_names{d}];
    print('-dpng',t_name), close
 
    Fig=figure(1);
    set(Fig,'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [40 20]);
    set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
    subplot(3,1,1), title('MP')
    pcolor(T,F,10*log10(Pw));shading flat, ylim([0 10]), caxis([-60 5])
    subplot(3,1,2), title('LF3')
    pcolor(T,F,10*log10(P3));shading flat, ylim([0 10]), caxis([-60 5])
    subplot(3,1,3), title('LF3 Prior')
    pcolor(T2,F2,10*log10(P32));shading flat, ylim([0 10]), caxis([-60 5])  
    t_name = ['liny_specgram3_compare_' f_names{d}];
    print('-dpng',t_name)
    subplot(3,1,1), set(gca,'yscale','log'), ylim([0.2 10])
    subplot(3,1,2), set(gca,'yscale','log'), ylim([0.2 10])
    subplot(3,1,3), set(gca,'yscale','log'), ylim([0.2 10])
    t_name = ['logy_specgram3_compare_' f_names{d}];
    print('-dpng',t_name), close
   
    clear lf*
end
