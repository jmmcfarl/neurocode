clear all
cd C:\WC_Germany\DUAL-MP_recordings
load dual_mp_dir
load desynch_detect\desynch_times
load UDS_dur_raw\UDS_raw_data

dsf = 8;
Fsd = 2016/dsf;
params.Fs = Fsd;
params.tapers = [4 7];
params.err = [2 .05];
params.fpass = [0 45];
window = [100 100];
niqf = 2016/2;
hcf1 = 100/niqf;
[b1,a1] = butter(2,hcf1,'low');

for d = 1:length(dir_array)
    
    disp(sprintf('session %d',d))
    cd(dir_array{d});
    pwd
    
    load used_data lf*
    
    %bandlimit signals
    down_8 = filtfilt(b1,a1,lf8);
        down_8 = downsample(down_8,dsf);
    down_8 = zscore(down_8);

     down_13 = filtfilt(b1,a1,lf13);
        down_13 = downsample(down_13,dsf);
    down_13 = zscore(down_13);
    
        down_15 = filtfilt(b1,a1,lf15);
        down_15 = downsample(down_15,dsf);
    down_15 = zscore(down_15);

        down_16 = filtfilt(b1,a1,lf16);
        down_16 = downsample(down_16,dsf);
    down_16 = zscore(down_16);

    if exist('lf14')
        down_14 = filtfilt(b1,a1,lf14);
        down_14 = downsample(down_14,dsf);
        down_14 = zscore(down_14);
    end
    
    desynch_start = desynch_start_times{d}*Fsd;
    desynch_stop = desynch_stop_times{d}*Fsd;
    
    if ~isempty(desynch_start)
        sMarkers = [[1;desynch_stop'] [desynch_start';length(down_8)]]
    else
        sMarkers = [1 length(down_8)];
    end
    
[C8_13(d,:),Phi,Smn,Smm,f,ConfC(d),PhiStd,C8_13err(d,:,:)] = coherencyc_unequal_length_trials([down_8 down_13],window, params, sMarkers);
 [C8_15(d,:),Phi,Smn,Smm,f,ConfC(d),PhiStd,C8_15err(d,:,:)] = coherencyc_unequal_length_trials([down_8 down_15],window, params, sMarkers);
[C8_16(d,:),Phi,Smn,Smm,f,ConfC(d),PhiStd,C8_16err(d,:,:)] = coherencyc_unequal_length_trials([down_8 down_16],window, params, sMarkers);
   
if exist('lf14')
    [C8_14(d,:),Phi,Smn,Smm,f,ConfC(d),PhiStd,C8_14err(d,:,:)] = coherencyc_unequal_length_trials([down_8 down_14],window, params, sMarkers);
end
%    
if exist('lf14')
        plot(f,C8_13(d,:),'linewidth',2)
    hold on
    plot(f,C8_14(d,:),'c','linewidth',2)
    plot(f,C8_15(d,:),'r','linewidth',2)
    plot(f,C8_16(d,:),'k','linewidth',2)
    line([0 45],[ConfC(d) ConfC(d)],'Color','k')
    xlim([0 45])
    legend('8-13','8-14','8-15','8-16')
    tname = ['C:\WC_Germany\DUAL-MP_recordings\coherence\wband_' f_names{d}];
    print('-dpng',tname);
    xlim([0 1])
    tname = ['C:\WC_Germany\DUAL-MP_recordings\coherence\8_' f_names{d}];
    print('-dpng',tname);
    close
else
    plot(f,C8_13(d,:),'linewidth',2)
    hold on
    plot(f,C8_15(d,:),'r','linewidth',2)
    plot(f,C8_16(d,:),'k','linewidth',2)
    line([0 45],[ConfC(d) ConfC(d)],'Color','k')
    xlim([0 45])
    legend('8-13','8-15','8-16')
    tname = ['C:\WC_Germany\DUAL-MP_recordings\coherence\wband_' f_names{d}];
    print('-dpng',tname);
    xlim([0 1])
    tname = ['C:\WC_Germany\DUAL-MP_recordings\coherence\8_' f_names{d}];
    print('-dpng',tname);
    close
end
   
    clear down* lf*
    
end


save C:\WC_Germany\DUAL-MP_recordings\coherence\coherence_data C8* f ConfC