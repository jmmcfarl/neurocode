clear all

load E:\WC_Germany\JMM_Analysis\file_directories_jmm.mat


%initialize overall variables
wlf2_cohere_overall = cell(1,length(dir_array));
wlf3_cohere_overall = cell(1,length(dir_array));
wlf5_cohere_overall = cell(1,length(dir_array));
wlf7_cohere_overall = cell(1,length(dir_array));
wlf8_cohere_overall = cell(1,length(dir_array));
wlf2_phase_overall = cell(1,length(dir_array));
wlf3_phase_overall = cell(1,length(dir_array));
wlf5_phase_overall = cell(1,length(dir_array));
wlf7_phase_overall = cell(1,length(dir_array));
wlf8_phase_overall = cell(1,length(dir_array));

d = 1;

while d <= length(dir_array)
    
    cd(dir_array{d})
    pwd

    
    load used_data_multLFP
    load used_data_multLFP_2
    
    
    Fs = median(CSC8_SampleFrequencies);
    params.Fs = Fs;
    params.err = [1 0.01];
    params.fpass = [0.01 100];
    win = 20; %window in seconds
    numSegs = 10;

    segLength = floor(length(synct)/numSegs);
    
    for n = 1:numSegs
        [C_L2(n,:),phi_L2(n,:),S12_L2(n,:),S1,S2,f,confC_L2(n),phistd_L2(n,:)] = coherencysegc(wcv_minus_spike((n-1)*segLength+1:n*segLength),lf2((n-1)*segLength+1:n*segLength),win,params);
        [C_L3(n,:),phi_L3(n,:),S12_L3(n,:),S1,S2,f,confC_L3(n),phistd_L3(n,:)] = coherencysegc(wcv_minus_spike((n-1)*segLength+1:n*segLength),lf3((n-1)*segLength+1:n*segLength),win,params);
        [C_L5(n,:),phi_L5(n,:),S12_L5(n,:),S1,S2,f,confC_L5(n),phistd_L5(n,:)] = coherencysegc(wcv_minus_spike((n-1)*segLength+1:n*segLength),lf5((n-1)*segLength+1:n*segLength),win,params);
        [C_L7(n,:),phi_L7(n,:),S12_L7(n,:),S1,S2,f,confC_L7(n),phistd_L7(n,:)] = coherencysegc(wcv_minus_spike((n-1)*segLength+1:n*segLength),lf7((n-1)*segLength+1:n*segLength),win,params);
        [C_L8(n,:),phi_L8(n,:),S12_L8(n,:),S1,S2,f,confC_L8(n),phistd_L8(n,:)] = coherencysegc(wcv_minus_spike((n-1)*segLength+1:n*segLength),lf8((n-1)*segLength+1:n*segLength),win,params);

    end
    
    %calculate overall coherency params
    C_overall_L2{d} = mean(C_L2,1);
    C_overall_L2_SE{d} = std(C_L2,[],1)/sqrt(numSegs);
    phi_overall_L2{d} = mean(phi_L2,1);
    phi_overall_L2_SE{d} = std(phi_L2,[],1)/sqrt(numSegs);
    S12_overall_L2{d} = mean(S12_L2,1);
    S12_overall_L2_SE{d} = std(S12_L2,[],1)/sqrt(numSegs);
    confC_overall_L2{d} = mean(confC_L2);
    phistd_overall_L2{d} = mean(phistd_L2,1);
    C_overall_L3{d} = mean(C_L3,1);
    C_overall_L3_SE{d} = std(C_L3,[],1)/sqrt(numSegs);
    phi_overall_L3{d} = mean(phi_L3,1);
    phi_overall_L3_SE{d} = std(phi_L3,[],1)/sqrt(numSegs);
    S12_overall_L3{d} = mean(S12_L3,1);
    S12_overall_L3_SE{d} = std(S12_L3,[],1)/sqrt(numSegs);
    confC_overall_L3{d} = mean(confC_L3);
    phistd_overall_L3{d} = mean(phistd_L3,1);
     C_overall_L5{d} = mean(C_L5,1);
    C_overall_L5_SE{d} = std(C_L5,[],1)/sqrt(numSegs);
    phi_overall_L5{d} = mean(phi_L5,1);
    phi_overall_L5_SE{d} = std(phi_L5,[],1)/sqrt(numSegs);
    S12_overall_L5{d} = mean(S12_L5,1);
    S12_overall_L5_SE{d} = std(S12_L5,[],1)/sqrt(numSegs);
    confC_overall_L5{d} = mean(confC_L5);
    phistd_overall_L5{d} = mean(phistd_L5,1);
    C_overall_L7{d} = mean(C_L7,1);
    C_overall_L7_SE{d} = std(C_L7,[],1)/sqrt(numSegs);
    phi_overall_L7{d} = mean(phi_L7,1);
    phi_overall_L7_SE{d} = std(phi_L7,[],1)/sqrt(numSegs);
    S12_overall_L7{d} = mean(S12_L7,1);
    S12_overall_L7_SE{d} = std(S12_L7,[],1)/sqrt(numSegs);
    confC_overall_L7{d} = mean(confC_L7);
    phistd_overall_L7{d} = mean(phistd_L7,1);
    C_overall_L8{d} = mean(C_L8,1);
    C_overall_L8_SE{d} = std(C_L8,[],1)/sqrt(numSegs);
    phi_overall_L8{d} = mean(phi_L8,1);
    phi_overall_L8_SE{d} = std(phi_L8,[],1)/sqrt(numSegs);
    S12_overall_L8{d} = mean(S12_L8,1);
    S12_overall_L8_SE{d} = std(S12_L8,[],1)/sqrt(numSegs);
    confC_overall_L8{d} = mean(confC_L8);
    phistd_overall_L8{d} = mean(phistd_L8,1);

    
    d = d +1

    save E:\WC_Germany\JMM_Analysis\overall_data_jmm_lfp_cohere.mat *overall* f d

    %prepare for next iteration
    clear all

    load E:\WC_Germany\JMM_Analysis\file_directories_jmm.mat
    load E:\WC_Germany\JMM_Analysis\overall_data_jmm_lfp_cohere.mat

    
end