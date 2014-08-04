clear all

load C:\WC_Germany\Persistent_activity\dir_tree_update
load C:\WC_Germany\Persistent_activity\desynch_extract\desynch_points

dsf = 8;
Fsd = 2016/dsf;
params.Fs = Fsd;
params.tapers = [4 7];
params.err = [2 .05];
niqf = 2016/2;
lcf = 0.05/niqf;
hcf = 40/niqf;
[b,a] = butter(2,[lcf hcf]);
maxlag = 10*Fsd;

for d = 1:length(dir_array)
    
    disp(sprintf('session %d',d))
    cd(dir_array{d});
    pwd
    
    load used_data wcv_minus_spike lf8
    
    %bandlimit signals
    down_w = filtfilt(b,a,wcv_minus_spike);
    down_8 = filtfilt(b,a,lf8);

    down_w = downsample(down_w,dsf);
    down_8 = downsample(down_8,dsf);
      
    %zscore
    down_w = zscore(down_w);
    down_8 = zscore(down_8);
    
    [w_acorr(d,:),lags] = xcov(down_w,maxlag,'coeff');
    [l_acorr(d,:),lags] = xcov(down_8,maxlag,'coeff');
    [x_cor(d,:),lags] = xcov(down_w,down_8,maxlag,'coeff');
    
        
    
    clear down_w down_8 wcv* lf8 
    
end


save C:\WC_Germany\Persistent_activity\corr_data w_acorr lags l_acorr Fsd x_cor