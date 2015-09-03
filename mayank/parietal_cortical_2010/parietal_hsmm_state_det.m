%inputs:
% raw_data
% raw_Fs
clear all
close all

addpath('G:\Code\smoothing\software')
addpath('G:\Code\FullBNT-1.0.4\KPMstats\')
addpath('G:\Code\FullBNT-1.0.4\netlab3.3')
addpath('G:\WC_Germany\new_stellate_analysis\')
addpath('G:\WC_Germany\hsmm_state_detection')
addpath('G:\WC_Germany\parietal_cortical_2010\')

cd G:\WC_Germany\parietal_cortical_2010\
load parietal_cortical_2010

raw_Fs = 2016;
dsf = 8;
niqf = raw_Fs/2;
Fsd = raw_Fs/dsf;
[b,a] = butter(2,[0.05/niqf 30/niqf]);
[b2,a2] = butter(2,[0.05/niqf 2/niqf]);

win_size = round(20*Fsd);

for d = 1:length(sess_data)
    
    cd(sess_data(d).directory)
    load used_data wcv_minus_spike lf8 lf5
    [hsmm_state_seqn,hmm_state_seq,hmm,Fs] = parietal_get_hsmm_state_seq(wcv_minus_spike,raw_Fs,sess_data(d).name,0);
    [hsmm_state_seq8n,hmm_state_seq8,hmm8,Fs] = parietal_get_hsmm_state_seq(lf8,raw_Fs,strcat(sess_data(d).name,'_lf8'),0);
%     save hsmm_state_seq hsmm* hmm* Fs
    
%     load hsmm_state_seq
%     
% hsmm_state_seq = hsmm_state_seqn;
% 
    wcv_f = filtfilt(b,a,wcv_minus_spike);
    wcv_f = zscore(downsample(wcv_f,dsf));
    wcv_f2 = filtfilt(b2,a2,wcv_minus_spike);
    wcv_f2 = zscore(downsample(wcv_f2,dsf));
    lf8_f = filtfilt(b,a,lf8);
    lf8_f = zscore(downsample(lf8_f,dsf));
    lf5_f = filtfilt(b,a,lf5);
    lf5_f = zscore(downsample(lf5_f,dsf));
    time = (1:length(wcv_f))/Fsd;
    stime = (1:length(hsmm_state_seqn))/Fs;
    nWins = round(length(time))/win_size;
%     subplot(2,1,1)
    plot(time,wcv_f), hold on
%     plot(time,wcv_f2,'r')
    plot(stime,hsmm_state_seqn,'k','linewidth',2)
%     plot(stime,hsmm_state_seq-2,'g','linewidth',2)
%     subplot(2,1,2)
    plot(time,lf8_f,'r'), hold on
    plot(time,lf5_f,'k')
    plot(stime,hsmm_state_seq8n-2,'g','linewidth',2)    
%     plot(stime,hsmm_state_seq8-2,'g','linewidth',2)    
    for w = 1:nWins
        begT = round((w-1)*win_size+1);
        endT= begT + win_size;
%         subplot(2,1,1)
        xlim([begT,endT]/Fsd)
        ylim([-3 3])
%         subplot(2,1,2)
%         xlim([begT,endT]/Fsd)
        pause
    end

    
    
end
