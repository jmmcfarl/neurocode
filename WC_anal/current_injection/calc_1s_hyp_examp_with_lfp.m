clear all
close all


load('C:\WC_Germany\current_injection\pyramidal\L3MEC_Pyramids_CIs_Pulsedata\2007-06-04_CWC_LFP_B');

Fs = 2e4;
dsf = 4;
Fsd = Fs/dsf;

mp_min = -1.1;
mp_max = 0.1;
nbins = 500;
mp_range = linspace(mp_min,mp_max,nbins);


send_times = find(time==1.9999);
sbeg_times = find(time==0);

num_sweeps = length(sbeg_times);
tvec = 0:1/Fsd:(2-1/Fsd);
sweep_mat = zeros(num_sweeps,length(tvec));
for i = 1:num_sweeps-1
    curdata = data(sbeg_times(i):send_times(i));
    sweep_mat(i,:) = downsample(curdata,dsf);
end

used_sweep = 24;
cur_sweep = sweep_mat(used_sweep,:);

datadir = 'C:\WC_Germany\EC_MPasci\A2007_06_04_CWC_LFP_B';
load(datadir)
data_string = datadir(25:end);
data_string = strcat(data_string,'_MP');
eval(['cur_data = ' data_string ';'])
prior_dist = gpkde(cur_data,.02,[mp_min;mp_max;nbins]);
mean_prior = mean(cur_data);
std_prior = std(cur_data);

cur_sweep = (cur_sweep - mean_prior)/std_prior;

plot(tvec,cur_sweep)
ylim([-3.5 4.5])
line([0.5 0.5],[-3.5 -3],'Color','k')
line([1.5 1.5],[-3.5 -3],'Color','k')
offset = 54.021 - 0.9495;
load('C:\WC_Germany\Entorhinal-WC_LFP\2007-06-04_CWC_LFP_B\2007-6-4_16-24-12\used_data')
mean_wcv = mean(wcv);
std_wcv = std(wcv);
load('C:\WC_Germany\Entorhinal-WC_LFP\2007-06-04_CWC_LFP_B\2007-6-4_16-47-27\used_data')
Fs = 2016;
niqf = Fs/2;
lcf = 0.05/niqf;
hcf = 4/niqf;
[b,a] = butter(2,[lcf hcf]);
lf8_d = filtfilt(b,a,lf8);

dsf = 2;
Fsd = Fs/dsf;
lf8_d = downsample(lf8_d,dsf);
lf8_d = zscore(lf8_d);
wcv = (wcv - mean_wcv)/std_wcv;
t = (1:length(wcv))/Fs;
t_d = (1:length(lf8_d))/Fsd;
figure
plot(t_d-offset,lf8_d,'r')
xlim([0 3])

figure
plot(t-offset,wcv)
xlim([0 3])
ylim([-3 5])



