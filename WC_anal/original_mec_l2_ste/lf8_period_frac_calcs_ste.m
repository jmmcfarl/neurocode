clear all

load C:\WC_Germany\JMM_analysis_ste\dir_tree_ste
load C:\WC_Germany\JMM_analysis_ste\UDS_dur_run_hist\data

dsf = 8;
Fsd = 2016/dsf;

for d = 1:length(dir_array)
    
   cd(dir_array{d})
   pwd
   
   load used_data lf8 
   
   
   lf8_d = downsample(lf8,dsf);
   
   lf8_period_f{d} = zeros(size(lf8_d));
      lf8_period_p{d} = zeros(size(lf8_d));

   for i = 1:length(up_trans8{d})-1
       
       period_samps1 = down_trans8{d}(i)-up_trans8{d}(i);
       period_samps2 = up_trans8{d}(i+1) - down_trans8{d}(i);
%       lf8_period_f{d}(up_trans8{d}(i)+1:down_trans8{d}(i+1)) = i+linspace(0,1,period_samps);
       lf8_period_f{d}(up_trans8{d}(i)+1:down_trans8{d}(i)) = i+linspace(1,period_samps1,period_samps1)/period_samps1*0.5;
       lf8_period_f{d}(down_trans8{d}(i)+1:up_trans8{d}(i+1)) = i+0.5+linspace(1,period_samps2,period_samps2)/period_samps2*0.5;
            lf8_period_p{d}(up_trans8{d}(i)+1:down_trans8{d}(i)) = linspace(1,period_samps1,period_samps1)/period_samps1*180;
       lf8_period_p{d}(down_trans8{d}(i)+1:up_trans8{d}(i+1)) = 180+linspace(1,period_samps2,period_samps2)/period_samps2*180;
 
   end
  
end

cd C:\WC_Germany\JMM_analysis_ste\
save lf8_period_f_data lf8_period_f lf8_period_p