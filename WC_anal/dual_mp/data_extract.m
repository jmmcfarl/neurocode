load all_eeg_data 

if unique(CSC8_SampleFrequencies) ~= 2016
    error('dropped packet')
end

[len,wid] = size(CSC8_Samples);

lf2 = reshape(CSC2_Samples,len*wid,1);
lf3 = reshape(CSC3_Samples,len*wid,1);
lf5 = reshape(CSC5_Samples,len*wid,1);
lf6 = reshape(CSC6_Samples,len*wid,1);
lf7 = reshape(CSC7_Samples,len*wid,1);
lf8 = reshape(CSC8_Samples,len*wid,1);
lf13 = reshape(CSC13_Samples,len*wid,1);
lf15 = reshape(CSC15_Samples,len*wid,1);
lf16 = reshape(CSC16_Samples,len*wid,1);

if exist('CSC14_Samples')
    lf14 = reshape(CSC14_Samples,len*wid,1);
    lf14 = detrend(lf14);
end
clear *_Samples

lf2 = detrend(lf2);
lf3 = detrend(lf3);
lf5 = detrend(lf5);
lf6 = detrend(lf6);
lf7 = detrend(lf7);
lf8 = detrend(lf8);
lf13 = detrend(lf13);
lf15 = detrend(lf15);
lf16 = detrend(lf16);

save used_data lf*
