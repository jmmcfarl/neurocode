clear all
close all

load C:\WC_Germany\JMM_Analysis_pyr\dir_tree_update
load pyr_heka_dir
load C:\WC_Germany\JMM_analysis_pyr\temp_sweep_times_d5

ot1 = 2;
ot2 = 2+3.2432;
sd1 = 49999/5e3;
sd2 = 33783/5e3;
h_offset = [6.3 6.3 5.2 1.8 4.4 4.5 4.3 3.6 3.6 4.5 3.8 6.5 5.8 5.7 5.6 3.7];
off_time = [ot1 ot1 ot1 ot1 ot1 ot2 ot1 ot2 ot2 ot2 ot2 ot2 ot2 ot2 ot2 ot2];
sweep_dur = [sd1 sd1 sd1 sd1 sd1 sd2 sd1 sd2 sd2 sd2 sd2 sd2 sd2 sd2 sd2 sd2];

d = 5;
delta = 1;
cd C:\WC_Germany\Persistent_activity\Final_figures\super_up_heka_traces

target_sweeps = find(sweep_times{d}(:,1) > 220 & sweep_times{d}(:,2) < 420);
time_axis = [];
data_vec = [];
for i = 1:length(target_sweeps)
    time_axis = [time_axis (sweep_times{d}(target_sweeps(i),1):2e-4:sweep_times{d}(target_sweeps(i),2))-222 nan(1,10e3)]; 
    data_vec = [data_vec A2007_05_30_CWC_LFP_C_MP(sweep_ids{d}(target_sweeps(i),1):sweep_ids{d}(target_sweeps(i),2))' nan(1,10e3)]; 
%     plot((sweep_times{d}(target_sweeps(i),1):2e-4:sweep_times{d}(target_sweeps(i),2))-220,A2007_05_30_CWC_LFP_C_MP(sweep_ids{d}(target_sweeps(i),1):...
%        sweep_ids{d}(target_sweeps(i),2)))
   
%    ylim([-0.7 0])
%    xlim([sweep_times{d}(target_sweeps(i),1)-delta-220 sweep_times{d}(target_sweeps(i),2)+delta-220])
%    print(sprintf('Trace%d',i),'-deps')
%    close
    
    
end