clear all

load C:\WC_Germany\Cortical_analysis\cortical_dir

dsf = 8;
Fsd = 2016/dsf;
niqf = 2016/2;
lcf = .05/niqf;
hcf = 40/niqf;
[b,a] = butter(2,[lcf hcf]);

for d = 1:length(over_dir)
    
    cd(over_dir{d})
    pwd
    
    load used_data lf8 lf7 lf6 lf5
    lf8_f = filtfilt(b,a,lf8);
    lf7_f = filtfilt(b,a,lf7);
    lf6_f = filtfilt(b,a,lf6);
    lf5_f = filtfilt(b,a,lf5);
    
    lf8_d = downsample(lf8_f,dsf);
    lf7_d = downsample(lf7_f,dsf);
    lf6_d = downsample(lf6_f,dsf);
    lf5_d = downsample(lf5_f,dsf);
    
    csd_u = lf8_d + lf6_d - 2*lf7_d;
    csd_l = lf7_d + lf5_d - 2*lf6_d;
    
    save csd_data csd*
    
    clear lf* csd*

end