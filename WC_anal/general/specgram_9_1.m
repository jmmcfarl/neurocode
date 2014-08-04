clear all

load E:\WC_Germany\JMM_Analysis\file_directories_jmm.mat

params.fpass = [0 10];
params.err = [2 .05];
params.tapers = [4 7];
params.Fs = 504;
window = [40 20];

for d = 1:length(dir_array)

    cd(dir_array{d})
    pwd


    load used_data lf8 wcv_minus_spike

    lf8_down = downsample(lf8,4);
    wcv_down = downsample(wcv_minus_spike,4);

    lf8_down_z = lf8_down - mean(lf8_down);
    lf8_down_z = lf8_down_z/std(lf8_down_z);

    wcv_down_z = wcv_down - mean(wcv_down);
    wcv_down_z = wcv_down_z/std(wcv_down_z);


    [S8,t,f,S8err] = mtspecgramc(lf8_down_z,window,params);
    [Sw,t,f,Swerr] = mtspecgramc(wcv_down_z,window,params);

    lf8_gr{d} = zeros(length(t),length(f));
    wcv_gr{d} = zeros(length(t),length(f));
    
    logS8 = 10*log10(S8);
    logSw = 10*log10(Sw);
    
    for i = 1:length(t)

        temp_freq = find(S8err(1,i,:)>Swerr(2,i,:));
        if ~isempty(temp_freq)
            lf8_gr{d}(i,temp_freq) = logS8(i,temp_freq)-logSw(i,temp_freq);
        end

        clear temp_freq

        temp_freq = find(Swerr(1,i,:)>S8err(2,i,:));
        if ~isempty(temp_freq)
            wcv_gr{d}(i,temp_freq) = logSw(i,temp_freq) - logS8(i,temp_freq);
        end
        clear temp_freq
        
    end
    
    Fig = figure(1)
    clf
    set(Fig,'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
    set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
    pcolor(t,f,lf8_gr{d}');shading flat;colorbar
    xlabel('Time (s)','FontSize',14)
    ylabel('Frequency (Hz)','FontSize',14)
    title('LF8 Greater than WCV','FontSize',16)
    ylim([0 2])
    tname = ['E:\WC_Germany\JMM_Analysis\specgram_91\lf8_gr_' f_names{d}];
    print('-dpng',tname);
    
    clf
    pcolor(t,f,wcv_gr{d}');shading flat;colorbar
    xlabel('Time (s)','FontSize',14)
    ylabel('Frequency (Hz)','FontSize',14)
    title('Wcv Greater than LF8','FontSize',16)
    ylim([0 2])
    tname = ['E:\WC_Germany\JMM_Analysis\specgram_91\wcv_gr_' f_names{d}];
    print('-dpng',tname);
    
    




end

