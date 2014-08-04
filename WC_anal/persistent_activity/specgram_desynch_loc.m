clear all

load C:\WC_Germany\Persistent_activity\dir_tree_update

dsf = 8;
Fsd = 2016/dsf;
niqf = 2016/2;
[b,a] = butter(2,0.1/niqf,'high');

for d = 1:17
    
    cd(dir_array{d})
    pwd
    
    load used_data wcv_minus_spike lf8
    
    wcv_f = filtfilt(b,a,wcv_minus_spike);
    lf8_f = filtfilt(b,a,lf8);
    wcv_f = downsample(wcv_f,dsf);
    lf8_f = downsample(lf8_f,dsf);
    wcv_f = zscore(wcv_f);
    lf8_f = zscore(lf8_f);
    
    
    
    [Sw,F,T] = SPECTROGRAM(wcv_f,round(Fsd*20),round(Fsd*2),[],Fsd);
    Sw = abs(Sw).^2;
    
    %find max power in UDS band
    uds_min_f = find(F > 0.2,1,'first');
    uds_max_f = find(F > 1,1,'first');
    max_uds_pow = max(Sw(uds_min_f:uds_max_f,:));
    
    %find max power in theta band
    theta_min_f = find(F > 2,1,'first');
    theta_max_f = find(F > 6,1,'first');
    max_theta_pow = max(Sw(theta_min_f:theta_max_f,:));
    
    [S8,F,T] = SPECTROGRAM(lf8_f,round(Fsd*20),round(Fsd*2),[],Fsd);
    S8 = abs(S8).^2;
    
     %find max power in UDS band
    max_uds_pow8 = max(S8(uds_min_f:uds_max_f,:));
    
    %find max power in theta band
    max_theta_pow8 = max(S8(theta_min_f:theta_max_f,:));

    
    figure
    subplot(2,1,1)
    pcolor(T,F,10*log10(Sw));shading flat
    ylim([0 5])
    caxis([0 60])
    subplot(2,1,2)
    plot(T,max_theta_pow./max_uds_pow)
    hold on
    plot(T,max_uds_pow,'r')
    
    figure
    subplot(2,1,1)
    pcolor(T,F,10*log10(S8));shading flat
     ylim([0 5])
     caxis([0 60])
     subplot(2,1,2)
     plot(T,max_theta_pow8./max_uds_pow8)
     hold on
     plot(T,max_uds_pow,'r')
     line([0 max(T)],[5e5 5e5],'Color','k')
     
     pause
     close all
   
    
end