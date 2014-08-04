clear all
close all
cd C:\WC_Germany\Cortical_analysis

load C:\WC_Germany\Cortical_analysis\cortical_dir
load C:\WC_Germany\Cortical_analysis\desynch_detect\desynch_times
load C:\WC_Germany\Cortical_analysis\UDS_dur_raw\UDS_raw_data
load C:\WC_Germany\Cortical_analysis\UDS_synch_state_dur\UDS_synch_state_dur_data

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

for d = 1:length(sess_data)
    
    disp(sprintf('session %d',d))
    cd(sess_data(d).directory);
    pwd
    
    load used_data wcv_minus_spike lf8 lf5
    
    %bandlimit signals
    down_w = filtfilt(b1,a1,wcv_minus_spike);
    down_8 = filtfilt(b1,a1,lf8);
      down_5 = filtfilt(b1,a1,lf5);
  
    down_w = downsample(down_w,dsf);
    down_8 = downsample(down_8,dsf);
    down_5 = downsample(down_5,dsf);
    
    %zscore
    down_w = zscore(down_w);
    down_8 = zscore(down_8);
    down_5 = zscore(down_5);
    
    desynch_start = desynch_start_times{d}*Fsd;
    desynch_stop = desynch_stop_times{d}*Fsd;
    
    if ~isempty(desynch_start)
        sMarkers = [[1;desynch_stop'] [desynch_start';length(down_w)]]
    else
        sMarkers = [1 length(down_w)];
    end
    
    slength = sMarkers(:,2) - sMarkers(:,1);
    sMarkers(slength < 500,:) = [];
    
[Cmn(d,:),Phimn(d,:),Smn,Smm,f,ConfC(d),PhiStd,Cerr(d,:,:)] = coherencyc_unequal_length_trials([down_w down_8],window, params, sMarkers);
[Cmn5(d,:),Phimn5(d,:),Smn,Smm,f,ConfC5(d),PhiStd,Cerr5(d,:,:)] = coherencyc_unequal_length_trials([down_w down_5],window, params, sMarkers);

%     
    plot(f,Cmn(d,:))
    hold on
    plot(f,Cmn5(d,:),'r')
    legend('MP-LF8','MP-LF5')
    line([0 45],[ConfC(d) ConfC(d)],'Color','k')
    xlim([0 45])
        ylim([0 1])
    cell_name = ['L' sess_data(d).layer '_' sess_data(d).cell_type '_' sess_data(d).region '_' sess_data(d).name];
    tname = ['C:\WC_Germany\Cortical_analysis\coherence\wband_' cell_name];
    print('-dpng',tname);
    close
 

    
    plot(f,Cmn(d,:),'linewidth',2)
    hold on
%     plot(f,squeeze(Cerr(d,1,:)),'--')
%     plot(f,squeeze(Cerr(d,2,:)),'--')
plot(f,Cmn5(d,:),'r','linewidth',2)
legend('MP-LF8','MP-LF5')
    line([0 1],[ConfC(d) ConfC(d)],'Color','k')
    xlim([0 1])
    cell_name = ['L' sess_data(d).layer '_' sess_data(d).cell_type '_' sess_data(d).region '_' sess_data(d).name];
    ylim([0 1])
    tname = ['C:\WC_Germany\Cortical_analysis\coherence\' cell_name];
    print('-dpng',tname);
    close
  
    clear down_w down_8 wcv* lf8 
    
end


save C:\WC_Germany\Cortical_analysis\coherence\coherence_data Cmn* f Phimn* Cerr* ConfC*