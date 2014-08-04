clear all

load C:\WC_Germany\JMM_Analysis_pyr\dir_tree_update
load C:\WC_Germany\JMM_Analysis_pyr\UDS_dur_run_hist_v2\UDS_dur_data
load C:\WC_Germany\JMM_analysis_pyr\spectra\specgram_data

for d = 1:length(dir_array)
    
    cur_max_C = max_val_C{d};
    cur_max_freq_8 = max_freq_S8{d};
    cur_max_freq_Sw = max_freq_Sw{d};
    
    cur_pow_Sw = max_val_Sw{d};
    cur_pow_S8 = max_val_S8{d};
    
    cur_max_up = max_up_dur{d};
    cur_max_up8 = max_up_dur8{d};
    
    cur_max_up(1:7) = [];
    cur_max_up(end-6:end) = [];
    cur_max_up8(1:7) = [];
    cur_max_up8(end-6:end) = [];

    [a,b] = corrcoef(cur_max_freq_8,cur_max_freq_Sw);
    c_max_f_8_max_f_w(d) = a(2,1);
    c_max_f_8_max_f_w_p(d) = b(2,1);

    [a,b] = corrcoef(cur_max_freq_Sw,cur_pow_S8);
    c_max_f_w_pow_8(d) = a(2,1);
    c_max_f_w_pow_8_p(d) = b(2,1);

    [a,b] = corrcoef(cur_pow_Sw,cur_pow_S8);
    c_pow_w_pow_8(d) = a(2,1);
    c_pow_w_pow_8_p(d) = b(2,1);

    
    [a,b] = corrcoef(cur_max_C,cur_max_freq_8);
    c_max_C_max_f_8(d) = a(2,1);
    c_max_C_max_f_8_p(d) = b(2,1);
    
    [a,b] = corrcoef(cur_max_C,cur_max_freq_Sw);
    c_max_C_max_f_w(d) = a(2,1);
    c_max_C_max_f_w_p(d) = b(2,1);
    
    [a,b] = corrcoef(cur_max_C,cur_pow_Sw);
    c_max_C_pow_w(d) = a(2,1);
    c_max_C_pow_w_p(d) = b(2,1);
    
    [a,b] = corrcoef(cur_max_C,cur_pow_S8);
    c_max_C_pow_8(d) = a(2,1);
    c_max_C_pow_8_p(d) = b(2,1);

    [a,b] = corrcoef(cur_max_C(~isnan(cur_max_up)),cur_max_up(~isnan(cur_max_up)));
    c_max_C_mup_w(d) = a(2,1);
    c_max_C_mup_w_p(d) = b(2,1);

    [a,b] = corrcoef(cur_max_freq_8(~isnan(cur_max_up)),cur_max_up(~isnan(cur_max_up)));
    c_max_f_8_mup_w(d) = a(2,1);
    c_max_f_8_mup_w_p(d) = b(2,1);
 
    [a,b] = corrcoef(cur_pow_S8(~isnan(cur_max_up)),cur_max_up(~isnan(cur_max_up)));
    c_pow_8_mup_w(d) = a(2,1);
    c_pow_8_mup_w_p(d) = b(2,1);
  
end