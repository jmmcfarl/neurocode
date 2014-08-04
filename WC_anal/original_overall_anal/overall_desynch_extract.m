%% locate periods of desynchronized activity and mark beginning and end
clear all

load('C:\WC_Germany\overall_calcs\overall_dir.mat')
load C:\WC_Germany\overall_calcs\UDS_dur_raw\UDS_raw_data_new_method

up_dur_thresh = 5; %length of maximum allowable lfp up state in s
remove_win = 15; %remove this much data on either side of offending up states

dsf = 8;
Fsd = 2016/dsf;

for d = 1:length(over_dir)
    
   desynch_start{d} = [];
   desynch_stop{d} = [];
   
   long_lfp_ups = find(up_state_dur8{d} > up_dur_thresh);
    
   if ~isempty(long_lfp_ups)
       
      temp_d_start = up_trans8{d}(long_lfp_ups) - remove_win*Fsd;
      temp_d_stop = down_trans8{d}(long_lfp_ups) + remove_win*Fsd;
       
      %make sure states dont overlap now
      bad_d_start = [];
      for i = 2:length(temp_d_start)
         
          if temp_d_start(i) < temp_d_stop(i-1)
              bad_d_start = [bad_d_start i];
          end
          
      end
      
      temp_d_start(bad_d_start) = [];
      temp_d_stop(bad_d_start-1) = [];
      
      desynch_start{d} = temp_d_start;
      desynch_stop{d} = temp_d_stop;
      
      
   end
    
   plot(up_trans8{d}/Fsd,up_state_dur8{d},'.')
   hold on
   if ~isempty(desynch_start{d})
       plot(desynch_start{d}/Fsd,ones(size(desynch_start{d}))*5,'k*')
       plot(desynch_stop{d}/Fsd,ones(size(desynch_stop{d}))*5,'r*')
   end
    t_names = ['C:\WC_Germany\overall_calcs\desynch_extract\new_method_' num2str(cell_type(d)) '_' over_names{d}]
   print('-dpng',t_names)
   close all
    
save C:\WC_Germany\overall_calcs\desynch_extract\desynch_points_new_method desynch*    
    d
end

