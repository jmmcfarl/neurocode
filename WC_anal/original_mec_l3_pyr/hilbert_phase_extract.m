clear all

load C:\WC_Germany\JMM_analysis_pyr\dir_tree_update
load C:\WC_Germany\JMM_analysis_pyr\UDS_dur_run_hist_v2\up_per_data_10_28

dsf = 8;
Fsd = 2016/dsf;
niqf = 2016/2;
lcf = 0.05/niqf;
hcf = 1/niqf;
hcf2 = 10/niqf;
[b,a] = butter(2,[lcf hcf]);
[b2,a2] = butter(2,[lcf hcf2]);

for d = 1:length(dir_array)
    
   cd(dir_array{d})
   pwd
   
   load used_data lf8 wcv_minus_spike
   
   lf8_f = filtfilt(b,a,lf8);
   wcv_f = filtfilt(b,a,wcv_minus_spike);
   
   lf8_d = downsample(lf8_f,dsf);
   wcv_d = downsample(wcv_f,dsf);
   
   lf8_f2 = filtfilt(b2,a2,lf8);
   wcv_f2 = filtfilt(b2,a2,wcv_minus_spike);
   
   lf8_d2 = downsample(lf8_f2,dsf);
   wcv_d2 = downsample(wcv_f2,dsf);
   
  x8 = hilbert(lf8_d);
  xw = hilbert(wcv_d);
  
  phase8 = acos(real(x8)./imag(x8));
  phasew = acos(real(xw)./imag(xw));
    
  t = (1:length(lf8_d))/Fsd;  
  
  
end