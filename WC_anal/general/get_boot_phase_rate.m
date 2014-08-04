function [boot_phase_rate] = get_boot_phase_rate(lf8_phase,time_sls,phase_hist,lags,avg_rate_tsu,num_samps)

net_rate = zeros(length(phase_hist),1);
net_samps = zeros(length(phase_hist),1);
phase_binwidth = phase_hist(2)-phase_hist(1);
for i = 1:num_samps
    
   %select a random phase bin
   cur_phase_bin = ceil(rand*length(phase_hist));
   
   %find all instances where lf8_phase is within this range
   if cur_phase_bin == 1
       phase_inst = find(lf8_phase > 0 & lf8_phase < phase_hist(2));
   elseif cur_phase_bin == length(phase_hist)
       phase_inst = find(lf8_phase > phase_hist(end-1) & lf8_phase <= 360);
   else
       phase_inst = find(lf8_phase > phase_hist(cur_phase_bin) - phase_binwidth/2 & lf8_phase < phase_hist(cur_phase_bin)+phase_binwidth/2);
   end
   
   %select a time at random from the phase instances
   if ~isempty(phase_inst)
       
      rand_inst = phase_inst(ceil(rand*length(phase_inst)));
      cur_time_inst = time_sls(rand_inst);
      cur_time_bin = find(lags > cur_time_inst,1,'first');
        if isempty(cur_time_bin)
            cur_time_bin = length(lags);
        end
       net_rate(cur_phase_bin) = net_rate(cur_phase_bin) + avg_rate_tsu(cur_time_bin);
       net_samps(cur_phase_bin) = net_samps(cur_phase_bin)+1;
   end
    
   
end


boot_phase_rate = net_rate./net_samps;



