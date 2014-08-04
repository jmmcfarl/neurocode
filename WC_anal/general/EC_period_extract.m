%% GET UP STATE AMPLITUDE BATCH

clear all

load E:\WC_Germany\JMM_Analysis\file_directories_jmm.mat


% %initialize overall variables
median_spectrum_wcv = zeros(43,1301);
specgram_wcv_vals = cell(1,43);
median_spectrum_lf8 = zeros(43,1301);
specgram_lf8_vals = cell(1,43);

d = 1;

while d <= length(dir_array)

    cd(dir_array{d})
    pwd

    load used_data CSC8_SampleFrequencies wcv_minus_spike dt synct lf8

    Fs = mean(CSC8_SampleFrequencies);

    %zscore
    wcv_z = wcv_minus_spike - mean(wcv_minus_spike);
    wcv_z = wcv_z/std(wcv_z);
    lf8 = lf8-mean(lf8);
    lf8 = lf8/std(lf8);

params.Fs = Fs;
params.fpass = [0 40];
params.err = 0;
movingwin = [20 5];

[S,t,f] = mtspecgramc(wcv_z,movingwin,params);
    
    median_spectrum_wcv(d,:) = median(S,1);
    specgram_wcv_vals{d} = S;
    
        [S,t,f] = mtspecgramc(lf8,movingwin,params);

        median_spectrum_lf8(d,:) = median(S,1);
        specgram_lf8_vals{d} = S;
        
    clear S 
    
    d = d +1

    save E:\WC_Germany\JMM_Analysis\overall_data_jmm_wcv_UDS_trans_period.mat f t d median* specgram*

    %prepare for next iteration
    clear all

    load E:\WC_Germany\JMM_Analysis\file_directories_jmm.mat
    load E:\WC_Germany\JMM_Analysis\overall_data_jmm_wcv_UDS_trans_period.mat f t d median* specgram*


end



% 
% 
% 
