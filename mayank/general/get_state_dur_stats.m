num_cells = length(up_trans)
for i = 1:num_cells
   mean_up_dur(i) = nanmean(up_state_dur{i});
   mean_up_dur8(i) = nanmean(up_state_dur8{i});
   mean_down_dur(i) = nanmean(down_state_dur{i});
   mean_down_dur8(i) = nanmean(down_state_dur8{i});
   
   median_up_dur(i) = nanmedian(up_state_dur{i});
   median_up_dur8(i) = nanmedian(up_state_dur8{i});
   median_down_dur(i) = nanmedian(down_state_dur{i});
   median_down_dur8(i) = nanmedian(down_state_dur8{i});
   
   std_up_dur(i) = nanstd(up_state_dur{i});
   std_up_dur8(i) = nanstd(up_state_dur8{i});
   std_down_dur(i) = nanstd(down_state_dur{i});
   std_down_dur8(i) = nanstd(down_state_dur8{i});
    
   all_dc = nan(length(up_trans{i}),1);
   for j = 1:length(up_state_dur{i})-1
        all_dc(j) = up_state_dur{i}(j)/(up_state_dur{i}(j)+down_state_dur{i}(j));
   end
   mean_dc(i) = nanmean(all_dc);
   median_dc(i) = nanmedian(all_dc);
   
     all_dc = nan(length(up_state_dur8{i}),1);
   for j = 1:length(up_trans8{i})-1
        all_dc(j) = up_state_dur8{i}(j)/(up_state_dur8{i}(j)+down_state_dur8{i}(j));
   end
    mean_dc8(i) = nanmean(all_dc);
   median_dc8(i) = nanmedian(all_dc);
 
end