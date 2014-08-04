function [lfp_data,lfp_timestamps,Fs] = load_lfp_data(filename,dsf)

load(filename);
raw_Fs = 1/FullV.samper;
Fs = raw_Fs/dsf;
hcf = 0.8*Fs/2;
[b,a] = butter(4,hcf/(raw_Fs/2),'low');

lfp_data = double(FullV.V)*FullV.intscale(1)/FullV.intscale(2);
lfp_data = filtfilt(b,a,lfp_data);
lfp_data = downsample(lfp_data,dsf);
n_pts = length(lfp_data);
lfp_timestamps = FullV.start:1/Fs:(FullV.start+(n_pts-1)/Fs);