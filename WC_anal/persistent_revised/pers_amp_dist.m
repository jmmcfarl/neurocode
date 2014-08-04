cd C:\WC_Germany\persistent_revised
clear all
close all

load pers_revised_dir

Fs=  2016;
niqf = Fs/2;
lcf = 0.05/niqf;
hcf = 40/niqf;
[b,a] = butter(2,[lcf hcf]);

amp_range = linspace(-4, 8, 600);
amp_range_f = linspace(-4, 4, 400);

for d = 1:28
    
    cd(dir_array{d})
    
    load used_data wcv wcv_minus_spike lf8
    
    wcv_f = filtfilt(b,a,wcv_minus_spike);
    lf8_f = filtfilt(b,a,lf8);
    wcv_f = zscore(wcv_f);
    wcvz = zscore(wcv);
    lf8_f = zscore(lf8_f);
    lf8 = zscore(lf8);
    
    amp_dist(d,:) = gpkde(wcvz, -3, [-4; 8; 600]);
    amp_dist_f(d,:) = gpkde(wcv_f, -3, [-4; 4; 400]);
    lf8_dist_f(d,:) = gpkde(lf8_f,-3,[-4; 4; 400]);
    lf8_dist(d,:) = gpkde(lf8,-3,[-4; 4; 400]);
    d
end

cd C:\WC_Germany\persistent_revised
save amp_dist amp_range* amp_dist* lf8_dist*