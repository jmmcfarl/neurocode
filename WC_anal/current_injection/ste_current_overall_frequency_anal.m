
clear all

load C:\WC_Germany\current_injection\stellate\ste_current_dir

niqf = 2016/2;
[b,a] = butter(2,[0.1/niqf 10/niqf]);
dsf = 8;

params.Fs = 2016/dsf;
params.tapers = [10 19];
params.err = [2 .05];
params.fpass = [0 20];

for d = 1:length(dir_array)
    cd(dir_array{d})
    pwd
    load used_data
    
    lf8_f = filtfilt(b,a,lf8);
    lf8_fd = zscore(downsample(lf8_f,dsf));
    lf3_f = filtfilt(b,a,lf3);
    lf3_fd = zscore(downsample(lf3_f,dsf));
    wcv_f = filtfilt(b,a,wcv_minus_spike);
    wcv_fd = zscore(downsample(wcv_f,dsf));
    
[Sw,f,Serrw]=mtspectrumc(wcv_fd,params);    
[S3,f,Serr3]=mtspectrumc(lf3_fd,params);    
[S8,f,Serr8]=mtspectrumc(lf8_fd,params);    

[C3,phi,S12,S1,S2,f,confC3,phistd,Cerr3]=coherencyc(wcv_fd,lf3_fd,params);
[C8,phi,S12,S1,S2,f,confC8,phistd,Cerr8]=coherencyc(wcv_fd,lf8_fd,params);

figure
    pcolor(t,f,C3');shading flat
    t_names = ['C:\WC_Germany\current_injection\stellate\cohgrams\lf3_' f_names{d}];
    print('-dpng',t_names)
    close all
    
                figure
    pcolor(t,f,C8');shading flat
    t_names = ['C:\WC_Germany\current_injection\stellate\cohgrams\lf8_' f_names{d}];
    print('-dpng',t_names)
    close all
    
end

