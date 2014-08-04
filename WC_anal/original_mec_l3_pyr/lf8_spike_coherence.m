clear all

load C:\WC_Germany\JMM_Analysis_pyr\dir_tree_update
load C:\WC_Germany\JMM_Analysis_pyr\UDS_dur_run_hist_v2\up_per_data_10_28
Fs = 2016;
dsf = 8;
Fsd = Fs/dsf;

params.Fs = Fsd;
params.fpass = [0 100];
params.err = [2 0.05];
win = 100;

for d = 1:length(dir_array)
    
   cd(dir_array{d});
   pwd
   
   load used_data lf8 
   load spike_time
   
   lf8 = downsample(lf8,dsf);
   
   spike_times = spkid/Fsd;
   
[C,phi,S12,S1,S2,f,zerosp,confC,phistd,Cerr]=coherencysegcpt(lf8,spike_times,win,params,1,1);    
    
        
end
