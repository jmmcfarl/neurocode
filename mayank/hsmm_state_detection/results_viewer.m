clear all
close all

cd G:\WC_Germany\persistent_revised\
load pers_revised_dir_10_29_09
load G:\WC_Germany\persistent_revised\hmmstatedet_dur_raw\hmmstate_det_data.mat

Fs = 2016;
dsf = 8;
niqf = Fs/2;
Fsd = Fs/dsf;
[b,a] = butter(2,[0.05/niqf 30/niqf]);
win_size = round(20*Fsd);

for d = 23:length(dir_array)
    
    temp_dir = dir_array{d};
    temp_dir(1) = 'G';
    cd(temp_dir)
    disp(temp_dir)
    load used_data wcv_minus_spike lf8
    load hsmm_state_seq_lf
    
    wcv_f = filtfilt(b,a,wcv_minus_spike);
    wcv_f = zscore(downsample(wcv_f,dsf));
    lf8_f = filtfilt(b,a,lf8);
    lf8_f = zscore(downsample(lf8_f,dsf));
    time = (1:length(lf8_f))/Fsd;
    stime = (1:length(hsmm_state_seqn))/Fs;
    
    state_seq = zeros(size(wcv_f));
    for i = 1:length(up_trans{d})
        state_seq(up_trans{d}(i):down_trans{d}(i)) = 1;
    end
    state_seq8 = zeros(size(lf8_f));
    for i = 1:length(up_trans8{d})
        state_seq8(up_trans8{d}(i):down_trans8{d}(i)) = 1;
    end
    
    
    nWins = round(length(time))/win_size;
    subplot(2,1,1)
    plot(time,wcv_f), hold on
    plot(stime,hsmm_state_seqn,'k','linewidth',2)
    plot(time,state_seq-1,'g','linewidth',2)
    subplot(2,1,2)
    plot(time,lf8_f), hold on
    plot(stime,hsmm_state_seq8n,'k','linewidth',2)
    plot(time,state_seq8-1,'g','linewidth',2)
    for w = 1:nWins
        begT = round((w-1)*win_size+1);
        endT= begT + win_size;
        subplot(2,1,1)
        xlim([begT,endT]/Fsd)
        subplot(2,1,2)
        xlim([begT,endT]/Fsd)
        pause
    end
    close
end