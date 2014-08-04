clear all
close all

load C:\WC_Germany\persistent_revised\pers_revised_dir
load C:\WC_Germany\persistent_revised\UDS_dur_raw\UDS_raw_data
load C:\WC_Germany\persistent_revised\UDS_synch_state_dur\UDS_synch_state_dur_data

Fs = 2016;
niqf = Fs/2;
lcf = 0.05/niqf;
hcf = 40/niqf;
[b,a] = butter(2,[lcf hcf]);
dsf = 8;
Fsd = Fs/dsf;

backlag = 5*Fsd;
forwardlag = 10*Fsd;
lags = (-backlag:forwardlag)/Fsd;

for d = 1:28
%  d=14   
    
    cd(dir_array{d})
    disp(num2str(d))
    
    load used_data lf8 wcv_minus_spike
    
    wcv_f = filtfilt(b,a,wcv_minus_spike);
    lf8_f = filtfilt(b,a,lf8);
    
    down_w = downsample(wcv_f,dsf);
    down_8 = downsample(lf8_f,dsf);
    
    down_w = zscore(down_w);
    down_8 = zscore(down_8);
    
    %initialize
    lf8_utrig_mp_mat = zeros(length(synch_ups8{d}),length(lags));
    lf8_utrig_lf8_mat = zeros(length(synch_ups8{d}),length(lags));
    lf8_dtrig_mp_mat = zeros(length(synch_downs8{d}),length(lags));
    lf8_dtrig_lf8_mat = zeros(length(synch_downs8{d}),length(lags));
   
    
    %calculate mp utrigs
    for i = 1:length(synch_ups8{d})
        
        if up_trans8{d}(synch_ups8{d}(i)) > backlag && ...
                length(down_w) - up_trans8{d}(synch_ups8{d}(i)) > forwardlag
           
            lf8_utrig_mp_mat(i,:) = down_w(up_trans8{d}(synch_ups8{d}(i))-backlag:...
                up_trans8{d}(synch_ups8{d}(i))+forwardlag);
            lf8_utrig_lf8_mat(i,:) = down_8(up_trans8{d}(synch_ups8{d}(i))-backlag:...
                up_trans8{d}(synch_ups8{d}(i))+forwardlag);

        else
            
            lf8_utrig_mp_mat(i,:) = nan;
            lf8_utrig_lf8_mat(i,:) = nan;
            
        end
        
    end
    
       %calculate mp dtrigs
    for i = 1:length(synch_downs8{d})
        
        if down_trans8{d}(synch_downs8{d}(i)) > backlag && ...
                length(down_w) - down_trans8{d}(synch_downs8{d}(i)) > forwardlag
           
            lf8_dtrig_mp_mat(i,:) = down_w(down_trans8{d}(synch_downs8{d}(i))-backlag:...
                down_trans8{d}(synch_downs8{d}(i))+forwardlag);
            lf8_dtrig_lf8_mat(i,:) = down_8(down_trans8{d}(synch_downs8{d}(i))-backlag:...
                down_trans8{d}(synch_downs8{d}(i))+forwardlag);

        else
            
            lf8_dtrig_mp_mat(i,:) = nan;
            lf8_dtrig_lf8_mat(i,:) = nan;
            
        end
        
    end
 
  
    lf8_utrig_mp(d,:) = nanmean(lf8_utrig_mp_mat);
    lf8_utrig_lf8(d,:) = nanmean(lf8_utrig_lf8_mat);
    lf8_dtrig_mp(d,:) = nanmean(lf8_dtrig_mp_mat);
    lf8_dtrig_lf8(d,:) = nanmean(lf8_dtrig_lf8_mat);
    

end

clear *mat
save C:\WC_Germany\persistent_revised\trig_avgs\lf8_trig_avg_data_wideband lags lf8*