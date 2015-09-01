function [boot_phase_rate,phase_time_dist] = get_boot_phase_rate_v2(lf8_phase,time_sls,phase_hist,lags,avg_rate_tsu,num_samps)

% net_rate = zeros(length(phase_hist),1);
% net_samps = zeros(length(phase_hist),1);
phase_binwidth = phase_hist(2)-phase_hist(1);
phase_time_dist = zeros(length(phase_hist),length(lags));

for i = 1:length(phase_hist)
    
   %find all instances where lf8_phase is in the current phase range
   if i == 1
       cur_phases_locs = find(lf8_phase < phase_hist(1) + phase_binwidth/2);
   elseif i == length(phase_hist)
       cur_phases_locs = find(lf8_phase > phase_hist(end)-phase_binwidth/2);
   else
       cur_phases_locs = find(lf8_phase > phase_hist(i)-phase_binwidth/2 & lf8_phase < phase_hist(i)+phase_binwidth/2);
   end
    
   %get the corresponding list of all times corresponding to these phase
   %instances
   if ~isempty(cur_phases_locs)
        tslu_list = time_sls(cur_phases_locs);
        phase_time_dist(i,:) = hist(tslu_list,lags);
        %take num_samp random elements of this list and find the
        %approximate average rate at this tslu
        cur_rate_list = zeros(num_samps,1);
        for n = 1:num_samps
            
            rand_time = tslu_list(ceil(rand*length(tslu_list)));
            time_bin = find(lags > rand_time,1,'first');
            if isempty(time_bin)
                time_bin = length(lags);
            end
            cur_rate_list(n) = avg_rate_tsu(time_bin);
            
        end
        
        boot_phase_rate(i) = nanmean(cur_rate_list);
   
   end
    
   
   
   
end


% boot_phase_rate = net_rate./net_samps;



