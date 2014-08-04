%analyze LFP state transitions using integrated high freq power

load C:\WC_Germany\JMM_Analysis_pyr\dir_tree_update

Fs = 2016;
niqf = 2016/2;
lcf = 20/niqf;
hcf = 100/niqf;
% binW = 5;

[b,a] = butter(2,[lcf hcf]);

for d = 1:length(dir_array)

    cd(dir_array{d})
    pwd

    load used_data lf8

    lf8_f = filtfilt(b,a,lf8);

    lf8_p = jmm_smooth_1d(lf8_f.^2,50);
    
%     num_bins = floor(length(lf8_f)/binW);

%     lf8_amp = zeros(num_bins,1);

%     for i = 1:num_bins
% 
%         beg_pt = (i-1)*binW+1;
%         end_pt = beg_pt+binW;
% 
%         lf8_amp(i) = std(lf8_f(beg_pt:end_pt));
% 
%     end
% 
%     lf8_p = jmm_smooth_1d(lf8_amp,30);
%     
%     lf8_d = downsample(lf8,binW);
    
end