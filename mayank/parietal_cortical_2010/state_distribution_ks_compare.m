
clear all
close all

addpath('F:\WC_Germany\parietal_cortical_2010\')
addpath('F:\WC_Germany\hsmm_state_detection\')

cd F:\WC_Germany\parietal_cortical_2010\
load parietal_cortical_2010
load F:\WC_Germany\parietal_cortical_2010\desynch_times_individual.mat


frontal = find_struct_field_vals(sess_data,'region','frontal');
prefrontal = find_struct_field_vals(sess_data,'region','prefrontal');
parietal = find_struct_field_vals(sess_data,'region','parietal');

% sess_data = sess_data(parietal);
% desynch_times = desynch_times(parietal);

%get rid of interneurons
interneurons = find_struct_field_vals(sess_data,'cell_type','interneuron');
sess_data(interneurons) = [];
% desynch_times(interneurons) = [];

raw_Fs = 2016;
dsf = 8;
Fsd = raw_Fs/dsf;

% thomas_el = find_struct_field_vals(sess_data,'thom_elec',1);
% sess_data = sess_data(thomas_el);
% desynch_times = desynch_times(thomas_el);

n = length(sess_data);

for d = 1:n
    
    cdir = sess_data(d).directory;
    cdir(1) = 'F';
    cd(cdir)
    pwd
    
    load hsmm_state_seq8_seg_lf_indivdesynch_dp
    [state_durations] = compute_state_durations_seg(hsmm_bbstate_seq8,Fsd);
    
    %for down state
    emp_pmf(d,1,:) = hist(state_durations{1},hmm8.dur_range);
    emp_pmf(d,1,end) = 0;
    emp_pmf(d,1,:) = emp_pmf(d,1,:)/sum(emp_pmf(d,1,:));
        
    %for up state    
    emp_pmf(d,2,:) = hist(state_durations{2},hmm8.dur_range);
    emp_pmf(d,2,end) = 0;
    emp_pmf(d,2,:) = emp_pmf(d,2,:)/sum(emp_pmf(d,2,:));
        
    load fixmean_state_seqs8
    [state_durations_f] = compute_state_durations_seg(fixmean_state_seq8,Fsd);
    
    %for down state
    empf_pmf(d,1,:) = hist(state_durations_f{1},hmm8.dur_range);
    empf_pmf(d,1,end) = 0;
    empf_pmf(d,1,:) = empf_pmf(d,1,:)/sum(empf_pmf(d,1,:));
    %for up state
    empf_pmf(d,2,:) = hist(state_durations_f{2},hmm8.dur_range);
    empf_pmf(d,2,end) = 0;
    empf_pmf(d,2,:) = empf_pmf(d,2,:)/sum(empf_pmf(d,2,:));
    
        
    load hsmm_state_seq8_seg_hf_indivdesynch_dp
    [state_durations_hf] = compute_state_durations_seg(hsmm_state_seq8_hf,hmm8_hf.Fs);
    
    %for down state    
    emphf_pmf(d,1,:) = hist(state_durations_hf{1},hmm8_hf.dur_range);
    emphf_pmf(d,1,end) = 0;
    emphf_pmf(d,1,:) = emphf_pmf(d,1,:)/sum(emphf_pmf(d,1,:));
    %for up state
    emphf_pmf(d,2,:) = hist(state_durations_hf{2},hmm8_hf.dur_range);
    emphf_pmf(d,2,end) = 0;
    emphf_pmf(d,2,:) = emphf_pmf(d,2,:)/sum(emphf_pmf(d,2,:));    
    
    [dummy,ks_8f(d,1)] = kstest2(state_durations{1},state_durations_f{1});
    [dummy,ks_8f(d,2)] = kstest2(state_durations{2},state_durations_f{2});
    [dummy,ks_8hf(d,1)] = kstest2(state_durations{1},state_durations_hf{1});
    [dummy,ks_8hf(d,2)] = kstest2(state_durations{2},state_durations_hf{2});
    
end

