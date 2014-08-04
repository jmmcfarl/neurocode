clear all
close all
cd C:\WC_Germany\current_injection\Simultaneous_LFP

part_rec{1} = 'C:\WC_Germany\april_09_data\2009-04-05_CWC_LFP\2009-4-5_18-26-31';
part_rec{2} = 'C:\WC_Germany\april_09_data\2009-04-07\2009-4-7-19';
part_rec{3} = 'C:\WC_Germany\april_09_data\2009-04-13_A\2009-4-13-18';
part_rec{4} = 'C:\WC_Germany\april_09_data\2009-04-21_CWC_LFP\2009-4-21_22-20-11';

part_name{1} = '2009-4-5';
part_name{2} = '2009-4-7';
part_name{3} = '2009-4-13';
part_name{4} = '2009-4-21';


Fs = 2016;
niqf = Fs/2;
dsf = 8;
Fsd = Fs/dsf;
lcf = 0.2/niqf;
hcf = 40/niqf;

[b,a] = butter(2,[lcf hcf]);

  cur_cell = 3;
  window_size = 30;
  
    cd(part_rec{cur_cell})
    pwd
    load raw_data
    wcv_d = filtfilt(b,a,wcv_minus_spike);
    wcv_d = zscore(downsample(wcv_d,dsf));
    lf8_d = filtfilt(b,a,lf8);
    lf8_d = zscore(downsample(lf8_d,dsf));
    lf3_d = filtfilt(b,a,lf3);
    lf3_d = zscore(downsample(lf3_d,dsf));
    
    t = (1:length(wcv_d))/Fsd;
    

    load Events_File
    load raw_sync_times
    load raw_spike_time_jmm
    spk_times = synct(spkid)/1e6;
    synct_d = downsample(synct,dsf)/1e6;
    rec_dur = synct_d(end);
        stim_on = Events_TimeStamps(Events_Nttls == -3)/1e6;
        dstim_on = [Inf diff(stim_on)];
        dstim_on2 = [dstim_on(2:end) Inf];
        pulse_start = stim_on(dstim_on > 0.2);
        pulse_stop = stim_on(dstim_on < 0.2 & dstim_on2 > 0.2);
       
          plot(synct_d,lf3_d);hold on
        cur_beg = 0;
        cur_end = cur_beg + window_size;
        while cur_beg < rec_dur
            cur_pulse_start = pulse_start(pulse_start > cur_beg & pulse_start < cur_end);
            cur_pulse_stop = pulse_stop(pulse_stop > cur_beg & pulse_stop < cur_end);
            xlim([cur_beg cur_end])
            for i = 1:length(cur_pulse_start)
                line([cur_pulse_start(i) cur_pulse_start(i)],[-2 2],'Color','g')
            end
            for i = 1:length(cur_pulse_stop)
                line([cur_pulse_stop(i) cur_pulse_stop(i)],[-2 2],'Color','r')
            end
            
            cur_spikes = spk_times(spk_times > cur_beg & spk_times < cur_end);
            plot(cur_spikes,ones(size(cur_spikes)),'k.')
            
            pause
            cur_beg = cur_beg + window_size;
            cur_end = cur_beg + window_size;
        end
        
    
    
    

