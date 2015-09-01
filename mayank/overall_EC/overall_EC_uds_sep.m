clear all

load F:\WC_Germany\overall_EC\overall_EC_dir.mat
addpath('F:\WC_Germany\hsmm_state_detection\')
addpath('F:\WC_Germany\parietal_cortical_2010\')

drive_letter = 'F';
dsf = 8;
raw_Fs = 2016;
Fsd = raw_Fs/dsf;

min_uds_time = 250;

for d = 1:length(sess_data)
    
    fprintf('session %d\n',d)
    cdir = sess_data(d).directory;
    cdir(1) = drive_letter;
    cd(cdir)
    
%     load used_data lf2 lf5 lf3
%     load desynch_times_lf2r
%     lf2 = lf2/sess_data(d).gains(2);
%     lf5 = lf5/sess_data(d).gains(5);
%     lf2_r = lf2-lf5;
%     
%     total_time = length(lf2)/raw_Fs;
%     desynch_time = sum(desynch_times_lf2r(:,2)-desynch_times_lf2r(:,1));
%     uds_time_2r(d) = total_time - desynch_time;
%     
%     if max(isnan(lf2_r)) > 0
%         bad_lf2r(d) = 1;
%     else
%         bad_lf2r(d) = 0;
%     end
    
%     if bad_lf2r(d) == 0 && uds_time_2r(d) > min_uds_time
%         
%         datalen = length(lf2_r);
%         load ec_hmm_state_seq2r
%         lf2r_state_seq = hmm_bbstate_seq2r;
%         
%         meandiff = [];
%         for i = 1:hmm2r.Nsegs
%             meandiff = [meandiff; hmm2r.state(2).meanfun{i}-hmm2r.state(1).meanfun{i}];
%         end
%         %         meandiff = mean(meandiff);
%         covar1 = hmm2r.state(1).var;
%         covar2 = hmm2r.state(2).var;
%         kl1 = gauss_kl_div(meandiff',covar1,covar2);
%         kl2 = gauss_kl_div(-meandiff',covar2,covar1);
%         lf2r_kl_lf(d) = kl1+kl2;
%     end
    
%     load desynch_times_lf3
%     if max(isnan(lf3)) > 0
%         bad_lf3(d) = 1;
%     else
%         bad_lf3(d) = 0;
%     end
%     total_time = length(lf3)/raw_Fs;
%     desynch_time = sum(desynch_times_lf3(:,2)-desynch_times_lf3(:,1));
%     uds_time_3(d) = total_time - desynch_time;
%     if bad_lf3(d) == 0 && uds_time_3(d) > min_uds_time
%         
%         datalen = length(lf3);
%         load ec_hmm_state_seq3
%         lf3_state_seq = hmm_bbstate_seq3;
%         
%         meandiff = [];
%         for i = 1:hmm3.Nsegs
%             meandiff = [meandiff; hmm3.state(2).meanfun{i}-hmm3.state(1).meanfun{i}];
%         end
%         %         meandiff = mean(meandiff);
%         covar1 = hmm3.state(1).var;
%         covar2 = hmm3.state(2).var;
%         kl1 = gauss_kl_div(meandiff',covar1,covar2);
%         kl2 = gauss_kl_div(-meandiff',covar2,covar1);
%         lf3_kl_lf(d) = kl1+kl2;
%     end
    
    
    load ./ec_hmm_state_seq
    mp_state_seq = hmm_bbstate_seq;
    meandiff = [];
    for i = 1:hmm.Nsegs
        meandiff = [meandiff; hmm.state(2).meanfun{i}-hmm.state(1).meanfun{i}];
    end
    uds_time_mp(d) = length(meandiff)/Fsd;
    %     meandiff = mean(meandiff);
    covar1 = hmm.state(1).var;
    covar2 = hmm.state(2).var;
    kl1 = gauss_kl_div(meandiff',covar1,covar2);
    kl2 = gauss_kl_div(-meandiff',covar2,covar1);
    mp_kl_lf(d) =  kl1+kl2;
    
    load ./ec_hmm_state_seq8
    lf8_state_seq = hmm_bbstate_seq8;
    meandiff = [];
    for i = 1:hmm8.Nsegs
        meandiff = [meandiff; hmm8.state(2).meanfun{i}-hmm8.state(1).meanfun{i}];
    end
    uds_time_lf8(d) = length(meandiff)/Fsd;
    %     meandiff = mean(meandiff);
    covar1 = hmm8.state(1).var;
    covar2 = hmm8.state(2).var;
    kl1 = gauss_kl_div(meandiff',covar1,covar2);
    kl2 = gauss_kl_div(-meandiff',covar2,covar1);
    lf8_kl_lf(d) =  kl1+kl2;
    
end

cd F:\WC_Germany\overall_EC\
save overall_EC_UDS_sep uds_* lf8_kl_lf mp_kl_lf

