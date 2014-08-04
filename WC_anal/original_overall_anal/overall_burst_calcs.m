  
clear all
close all

load('C:\WC_Germany\overall_calcs\overall_dir.mat')

Fs = 2016;

isi_range = linspace(0,100);

for d = 1:length(over_dir)
    
    cd(over_dir{d})
    disp(num2str(d))

if exist('spike_time.mat')
        load spike_time   
    else 
        load spike_time_br
end

num_spikes(d) = length(spkid);

    isis = [0;diff(spkid)]/Fs*1000;
    isi_hist(d,:) = hist(isis,isi_range);
    isi_hist(d,:) = isi_hist(d,:)/sum(isi_hist(d,:));
    
        
      bar(isi_range,isi_hist(d,:))
      xlim([0 99]);
      title(['total spikes: ' num2str(num_spikes(d))])
      t_names = ['C:\WC_Germany\overall_calcs\burst\' num2str(cell_type(d)) '_' over_names{d}];
        print('-dpng',t_names);
        close
        
      
      
end

save C:\WC_Germany\overall_calcs\burst\burst_data isi_hist num_spikes isi_range