clear all
cd C:\WC_Germany\DUAL-MP_recordings
load dual_mp_dir
load desynch_detect\desynch_times

dsf = 8;
Fsd = 2016/dsf;
params.Fs = Fsd;
params.tapers = [4 7];
params.err = [2 .05];
params.fpass = [0 45];
window = [100 100];
niqf = 2016/2;
lcf = 0.1/niqf;
hcf1 = 100/niqf;
[b1,a1] = butter(2,[lcf hcf1]);

for d = 1:length(dir_array)
    
    disp(sprintf('session %d',d))
    cd(dir_array{d});
    pwd
    
    load used_data lf*
    
    %bandlimit signals
    down_8 = filtfilt(b1,a1,lf8);
    down_8 = downsample(down_8,dsf);
%     down_8 = zscore(down_8);
    
    down_13 = filtfilt(b1,a1,lf13);
    down_13 = downsample(down_13,dsf);
%     down_13 = zscore(down_13);
    
    down_15 = filtfilt(b1,a1,lf15);
    down_15 = downsample(down_15,dsf);
%     down_15 = zscore(down_15);
    
    down_16 = filtfilt(b1,a1,lf16);
    down_16 = downsample(down_16,dsf);
%     down_16 = zscore(down_16);
  
    if exist('lf14')
            down_14 = filtfilt(b1,a1,lf14);
    down_14 = downsample(down_14,dsf);
%     down_14 = zscore(down_14);
    end
    
    desynch_start = round(desynch_start_times{d}*Fsd);
    desynch_stop = round(desynch_stop_times{d}*Fsd);
    
    if ~isempty(desynch_start)
        sMarkers = [[1;desynch_stop'] [desynch_start';length(down_8)]]
    else
        sMarkers = [1 length(down_8)];
    end
    
    [ S8(d,:), f, S8err(d,:,:) ]= mtspectrumc_unequal_length_trials(down_8,window, params, sMarkers);
     [ S13(d,:), f, S13err(d,:,:) ]= mtspectrumc_unequal_length_trials(down_13,window, params, sMarkers);
    [ S15(d,:), f, S15err(d,:,:) ]= mtspectrumc_unequal_length_trials(down_15,window, params, sMarkers);
    [ S16(d,:), f, S16err(d,:,:) ]= mtspectrumc_unequal_length_trials(down_16,window, params, sMarkers);
 if exist('lf14')
         [ S14(d,:), f, S14err(d,:,:) ]= mtspectrumc_unequal_length_trials(down_14,window, params, sMarkers);
 end
  
 if exist('lf14')
        Fig = figure(1)
    clf
    set(Fig,'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
    set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
    plot(f,S8(d,:),'Color',[0.8 0.8 0.7])
    hold on
%     plot(f,10*log10(S13(d,:)))
%     plot(f,10*log10(S14(d,:)),'c')
%     plot(f,10*log10(S15(d,:)),'r')
%     plot(f,10*log10(S16(d,:)),'k')
    plot(f,S13(d,:))
    plot(f,S14(d,:),'c')
    plot(f,S15(d,:),'r')
    plot(f,S16(d,:),'k')

    xlim([0.1 10])
%     ylim([-20 5])
    legend('LF8','LF13','LF14','LF15','LF16')
    set(gca,'xscale','log')
        tname = ['C:\WC_Germany\DUAL-MP_recordings\power_spectra\log_' f_names{d}];
    print('-dpng',tname);
close
 
 else
    Fig = figure(1)
    clf
    set(Fig,'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
    set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
%     plot(f,10*log10(S8(d,:)),'Color',[0.8 0.8 0.7])
%     hold on
%     plot(f,10*log10(S13(d,:)))
%     plot(f,10*log10(S15(d,:)),'r')
%     plot(f,10*log10(S16(d,:)),'k')
    plot(f,S8(d,:),'Color',[0.8 0.8 0.7])
    hold on
    plot(f,S13(d,:))
    plot(f,S15(d,:),'r')
    plot(f,S16(d,:),'k')

    xlim([0.1 10])
%         ylim([-20 5])
    legend('LF8','LF13','LF15','LF16')
    set(gca,'xscale','log')
        tname = ['C:\WC_Germany\DUAL-MP_recordings\power_spectra\log_' f_names{d}];
    print('-dpng',tname);
close
 end 
     clear down* lf*

end


save C:\WC_Germany\DUAL-MP_recordings\power_spec_data S* f 