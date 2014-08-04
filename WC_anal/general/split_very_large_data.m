% Split large data files into smaller chunks

np = 2;

for k = 1:np
    load all_eeg_data CSC1_Samples CSC1_TimeStamps CSC8_Samples CSC8_TimeStamps 
    [len,wid1] = size(CSC1_Samples);
    mid1 = round(wid1/np);
    CSC1_Samples = CSC1_Samples(1:len,(k-1)*mid1+1:(k-1)*mid1+mid1-1);
    CSC1_TimeStamps = CSC1_TimeStamps((k-1)*mid1+1:(k-1)*mid1+mid1-1);
    
%     [len,wid1] = size(CSC8_Samples);
%     mid1 = round(wid1/np);
%     CSC8_Samples = CSC8_Samples(1:len,(k-1)*mid1+1:(k-1)*mid1+mid1-1);
%     CSC8_TimeStamps = CSC8_TimeStamps((k-1)*mid1+1:(k-1)*mid1+mid1-1);
    
    
    load all_eeg_data CSC2_Samples CSC3_Samples CSC4_Samples ...
        CSC5_Samples CSC6_Samples CSC7_Samples CSC8_Samples CSC8_TimeStamps ...
        CSC1_SampleFrequencies CSC8_SampleFrequencies CSC2_TimeStamps

% load all_eeg_data CSC1_SampleFrequencies CSC8_SampleFrequencies

    if k == 1
        !mkdir part1
        cd part1; save part1_eeg_data
        cd ..;
    elseif k == 2
        !mkdir part2
        cd part2; save part2_eeg_data
        cd ..;
    elseif k == 3;
        cd part3; save part3_eeg_data
        cd ..
    elseif k == 4;
        cd part3; save part4_eeg_data
    end
    clear CSC1_Samples CSC1_TimeStamps CSC2_Samples CSC3_Samples CSC4_Samples ...
        CSC5_Samples CSC6_Samples CSC7_Samples CSC8_Samples CSC8_TimeStamps ...
        CSC1_SampleFrequencies CSC8_SampleFrequencies CSC2_TimeStamps
end
