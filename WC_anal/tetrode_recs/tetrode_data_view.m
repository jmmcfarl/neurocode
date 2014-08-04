close all
window_size = 20;

num_cells = length(TT_times);
cmap = [1 0 0; 1 0.5 0; 0.5 0.5 0; 0 0.5 1; 0.5 0.5 0.5; 1 0 0.5];
Fs = 2016;
niqf = Fs/2;
dsf = 4;
Fsd = Fs/dsf;

lcf = 0.1/niqf;
hcf = 2/niqf;

[b,a] = butter(2,[lcf hcf]);

lf1_f = filtfilt(b,a,lf1_data);
lf12_f = filtfilt(b,a,lf12_data);

lf1_d = downsample(lf1_f,dsf);
lf12_d = downsample(lf12_f,dsf);
lf1_d = zscore(lf1_d);
lf12_d = zscore(lf12_d);

t = (1:length(lf1_d))/Fsd;
plot(t,lf1_d,'linewidth',2)
hold on
plot(t,lf12_d,'r','linewidth',2)
legend('MEC LFP','Cortical LFP')
for i = 1:num_cells
%     
    cur_spikes{i} = TT_times{i};
    cur_spikes{i} = cur_spikes{i}-time_stamps(1)+1/Fsd;
    cur_spikes{i} = cur_spikes{i}/1e6; %convert to s
    
    plot(cur_spikes{i},ones(size(cur_spikes{i}))*(0.5-0.2*(i-1)),'.','Color',cmap(i,:),'MarkerSize',6)
    mean_rate(i) = length(cur_spikes{i})/max(t);
%     
end

spike_hist = hist(cur_spikes{num_cells},t);
spike_hist_sm = jmm_smooth_1d(spike_hist,50)*4;

plot(t,spike_hist_sm,'k','linewidth',2)

mean_rate
%     
total_dur = max(t);
% 
num_wins = ceil(total_dur/window_size);

for i = 1:num_wins
   
    xlim([(i-1)*window_size i*window_size])
    ylim([-3 6])
    pause
    
end
    
    