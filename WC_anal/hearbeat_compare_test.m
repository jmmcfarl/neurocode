cd C:\WC_Germany\april_09_data\2009-04-13_A\2009-4-13-18
load 2009-4-13_18-51-42_all_eeg_data
addpath('C:\Documents and Settings\Administrator\Desktop\Mehtalab_Analysis\Zach')
Fs = median(CSC1_SampleFrequencies);
wcv = reshape(CSC1_Samples,1,numel(CSC1_Samples));
wcv = -wcv;
clear CSC1_Samples
wcv_timestamps = interpolateTimestamps(CSC1_TimeStamps, 512, length(wcv));

load ../2009-4-13_18-51-42_events
hbeat = Events_TimeStamps(find(Events_Nttls==-2));

dsf = 20;
Fsd = Fs/dsf;

wcv_d = downsample(wcv,dsf);
wcv_td = downsample(wcv_timestamps,dsf);
