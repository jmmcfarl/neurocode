clear all

load C:\WC_Germany\current_injection\stellate\ste_current_dir

niqf = 2016/2;
[b,a] = butter(2,[0.1/niqf 10/niqf]);

params.Fs = 2016/4;
params.tapers = [5 9];
params.err = [2 .05];
params.fpass = [0 10];
movingwindow = [25 2];

for d = 1:length(dir_array)
    cd(dir_array{d})
    pwd
    load used_data
    
    lf8_f = filtfilt(b,a,lf8);
    lf8_fd = zscore(downsample(lf8_f,4));
    lf3_f = filtfilt(b,a,lf3);
    lf3_fd = zscore(downsample(lf3_f,4));
    wcv_f = filtfilt(b,a,wcv_minus_spike);
    wcv_fd = zscore(downsample(wcv_f,4));
    
    [C3,phi,S12,S1,S2,t,f,confC3,phistd,Cerr]=cohgramc(wcv_fd,lf3_fd,[25 2],params);
    [C8,phi,S12,S1,S2,t,f,confC8,phistd,Cerr]=cohgramc(wcv_fd,lf8_fd,[25 2],params);
    
    C3(C3 < confC3) = nan;
    C8(C8 < confC8) = nan;
    
%     [S3,t,f,Serr3]=mtspecgramc(lf3_fd,[25 2],params);
%     [S8,t,f,Serr8]=mtspecgramc(lf8_fd,[25 2],params);
%     [Sw,t,f,Serrw]=mtspecgramc(wcv_fd,[25 2],params);
    
%     figure
%     pcolor(t,f,10*log10(S8'));shading flat
%     t_names = ['C:\WC_Germany\current_injection\stellate\specgrams\lf8_' f_names{d}];
%     print('-dpng',t_names)
%     close all
%    
%     figure
%     pcolor(t,f,10*log10(S3'));shading flat
%     t_names = ['C:\WC_Germany\current_injection\stellate\specgrams\lf3_' f_names{d}];
%     print('-dpng',t_names)
%     close all
% 
%         figure
%     pcolor(t,f,10*log10(Sw'));shading flat
%     t_names = ['C:\WC_Germany\current_injection\stellate\specgrams\wcv_' f_names{d}];
%     print('-dpng',t_names)
%     close all

            figure
    pcolor(t,f,C3');shading flat
    t_names = ['C:\WC_Germany\current_injection\stellate\cohgrams\lf3_v2_' f_names{d}];
    print('-dpng',t_names)
    close all
    
                figure
    pcolor(t,f,C8');shading flat
    t_names = ['C:\WC_Germany\current_injection\stellate\cohgrams\lf8_v2' f_names{d}];
    print('-dpng',t_names)
    close all
    
end