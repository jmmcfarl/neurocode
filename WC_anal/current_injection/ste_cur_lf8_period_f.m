clear all

load C:\WC_Germany\current_injection\stellate\ste_current_dir
load C:\WC_Germany\current_injection\stellate\UDS_dur_raw\UDS_raw_data

dsf = 1;
Fsd = 2016/dsf;

for d = 1:length(dir_array)
    
   cd(dir_array{d})
   pwd
   
   load used_data lf8 
   
   up_trans8{d} = up_trans8{d}*8;
   down_trans8{d} = down_trans8{d}*8;
   
   lf8_d = downsample(lf8,dsf);
   
   lf8_up_fract = 0.4
   lf8_down_phase = 144
   
   lf8_period_f{d} = zeros(size(lf8_d));
   lf8_period_p{d} = zeros(size(lf8_d));
   for i = 1:length(up_trans8{d})-1
       
       period_samps1 = down_trans8{d}(i)-up_trans8{d}(i);
       period_samps2 = up_trans8{d}(i+1) - down_trans8{d}(i);
       
       lf8_period_f{d}(up_trans8{d}(i)+1:down_trans8{d}(i)) = ...
           i+linspace(1,period_samps1,period_samps1)/period_samps1*lf8_up_fract;
       lf8_period_f{d}(down_trans8{d}(i)+1:up_trans8{d}(i+1)) = ...
           i+lf8_up_fract+linspace(1,period_samps2,period_samps2)/period_samps2*(1-lf8_up_fract);
       
       lf8_period_p{d}(up_trans8{d}(i)+1:down_trans8{d}(i)) = ...
           linspace(1,period_samps1,period_samps1)/period_samps1*lf8_down_phase;
       lf8_period_p{d}(down_trans8{d}(i)+1:up_trans8{d}(i+1)) = ...
           lf8_down_phase+linspace(1,period_samps2,period_samps2)/period_samps2*(360-lf8_down_phase);
   
   end
  
end

cd C:\WC_Germany\current_injection\stellate\
save lf8_period_f lf8_period_f lf8_period_p