clear all

load C:\WC_Germany\JMM_analysis_pyr\dir_tree_update
load C:\WC_Germany\JMM_analysis_pyr\UDS_dur_run_hist_v2\up_per_data_10_28

dsf = 8;
Fsd = 2016/dsf;

for d = 1:length(dir_array)
    
   cd(dir_array{d})
   pwd
   
   load used_data lf8 
   
   
   lf8_d = downsample(lf8,dsf);
   
   lf8_period_f{d} = zeros(size(lf8_d));
   
   for i = 1:length(up_trans8{d})-1
       
       period_samps = up_trans8{d}(i+1)-up_trans8{d}(i);
       
      lf8_period_f{d}(up_trans8{d}(i)+1:up_trans8{d}(i+1)) = i+linspace(0,1,period_samps);
       
   end
  
end