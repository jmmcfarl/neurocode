%% for pyramidal cells
clear all
load C:\WC_Germany\JMM_Analysis_pyr\dir_tree_update
load pyr_heka_dir

ot1 = 2;
ot2 = 2+3.2432;
sd1 = 49999/5e3;
sd2 = 33783/5e3;

h_offset = [6.3 6.3 5.2 1.8 4.4 4.5 4.3 3.6 3.6 4.5 3.8 6.5 5.8 5.7 5.6 3.7];
off_time = [ot1 ot1 ot1 ot1 ot1 ot2 ot1 ot2 ot2 ot2 ot2 ot2 ot2 ot2 ot2 ot2];
sweep_dur = [sd1 sd1 sd1 sd1 sd1 sd2 sd1 sd2 sd2 sd2 sd2 sd2 sd2 sd2 sd2 sd2];

% for d = 1:length(f_loc)
    d=5
    cd(dir_array{d})
    pwd
    
    load(f_loc{d})
    dat_name = [f_loc{d} '_MP'];
    dat_name(1:24) = [];
    eval(['detrend(' dat_name ');'])       
      
    load used_data wcv
    
   t_axis = [1:length(wcv)]/2016;
   eval(['h_length = length(' dat_name ');'])
   h_secs = h_length/5e3;
   num_sweeps = floor(h_secs/10);
   
   %cycle through sweeps and add 
   sweep_times{d} = zeros(num_sweeps+1,2);
   sweep_ids{d} = zeros(num_sweeps+1,2);
    for n = 1:num_sweeps
        sweep_times{d}(n,1) = h_offset(d)+(n-1)*sweep_dur(d)+(n-1)*off_time(d);
        sweep_times{d}(n,2) = sweep_times{d}(n,1)+sweep_dur(d);
        sweep_ids{d}(n,1) = round((n-1)*sweep_dur(d)*5e3+1);
        sweep_ids{d}(n,2) = round(sweep_ids{d}(n,1)+sweep_dur(d)*5e3);
    end
    
    sweep_ids{d}(n+1,1) = sweep_ids{d}(n,2) + 1;
    sweep_ids{d}(n+1,2) = h_length;
    extra_sig = h_length - sweep_ids{d}(n+1,1);
    sweep_times{d}(n+1,1) = sweep_times{d}(n,2)+off_time(d);
    sweep_times{d}(n+1,2) = sweep_times{d}(n+1,1)+extra_sig/5e3;
      
% end
% 
% save C:\WC_Germany\JMM_analysis_pyr\heka_anal\sweep_times