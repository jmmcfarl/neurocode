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

cd G:\WC_Germany\persistent_revised\
load pers_revised_dir_10_29_09
raw_Fs = 2016;
dsf = 8;
niqf = raw_Fs/2;
Fsd = raw_Fs/dsf;
[b,a] = butter(2,[0.05/niqf 30/niqf]);
[b2,a2] = butter(2,[0.05/niqf 2/niqf]);

win_size = round(20*Fsd);

for d = 1:length(dir_array)
    
    temp_dir = dir_array{d};
    temp_dir(1) = 'G';
    cd(temp_dir)
    load used_data wcv_minus_spike lf8 lf7
    [hsmm_state_seq,hmm_state_seq,hmm,Fs] = get_hsmm_state_seq(wcv_minus_spike,raw_Fs,f_names{d});
%     [hsmm_state_seq8,hmm_state_seq8,hmm8,Fs] = get_hsmm_state_seq(lf8,raw_Fs,strcat(f_names{d},'_lf8'));
%     save hsmm_state_seq_hf hsmm* hmm* Fs
    
    load hsmm_state_seq
    
    wcv_f = filtfilt(b,a,wcv_minus_spike);
    wcv_f = zscore(downsample(wcv_f,dsf));
%     wcv_f2 = filtfilt(b2,a2,wcv_minus_spike);
%     wcv_f2 = zscore(downsample(wcv_f2,dsf));
    lf8_f = filtfilt(b,a,lf8);
    lf8_f = zscore(downsample(lf8_f,dsf));
     lf7_f = filtfilt(b,a,lf7);
    lf7_f = zscore(downsample(lf7_f,dsf));
   time = (1:length(lf8_f))/Fsd;
    stime = (1:length(hsmm_state_seq))/Fs;
    nWins = round(length(time))/win_size;
%     subplot(2,1,1)
%     plot(time,wcv_f), hold on
%     plot(time,wcv_f2,'r')
%     plot(stime,hsmm_state_seq,'k','linewidth',2)
%     plot(stime,hsmm_state_seq-2,'g','linewidth',2)
%     subplot(2,1,2)
    plot(time,lf8_f), hold on
    plot(time,lf5_f,'r')
    plot(stime,hsmm_state_seq8,'k','linewidth',2)    
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
% 
    
    
end
