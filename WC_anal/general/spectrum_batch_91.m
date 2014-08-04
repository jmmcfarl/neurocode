clear all

load E:\WC_Germany\JMM_Analysis\file_directories_jmm.mat

params.fpass = [0 10];
params.err = [2 .05];
params.tapers = [10 17];
dsf = 8;
params.Fs = 2016/dsf;

for d = 1:36

    cd(dir_array{d})
    pwd


    load used_data lf8 wcv_minus_spike

    lf8_down = downsample(lf8,dsf);
    wcv_down = downsample(wcv_minus_spike,dsf);

    lf8_down_z = lf8_down - mean(lf8_down);
    lf8_down_z = lf8_down_z/std(lf8_down_z);

    wcv_down_z = wcv_down - mean(wcv_down);
    wcv_down_z = wcv_down_z/std(wcv_down_z);


    [S8{d},f{d},varS,C,S8err{d}]=mtspectrumsegc(lf8_down_z,[100 10],params);
    [Sw{d},f{d},varS,C,Swerr{d}]=mtspectrumsegc(wcv_down_z,[100 10],params);
        
    Fig = figure(1)
    clf
    set(Fig,'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
    set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
    plot(f{d},10*log10(S8{d}),f{d},10*log10(S8err{d}(1,:)),'b--',f{d},10*log10(S8err{d}(2,:)),'b--')
    hold on
    plot(f{d},10*log10(Sw{d}),'r',f{d},10*log10(Swerr{d}(1,:)),'r--',f{d},10*log10(Swerr{d}(2,:)),'r--')
    legend('Lf8','Wcv')
    xlabel('Time (s)','FontSize',14)
    ylabel('Frequency (Hz)','FontSize',14)
    title('LF8 Greater than WCV','FontSize',16)
    xlim([0 2])
    tname = ['E:\WC_Germany\JMM_Analysis\spectrum_total\tot_spectrum_' f_names{d}];
    print('-dpng',tname);
    
    
    




end

