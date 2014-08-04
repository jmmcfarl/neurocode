function [jit_avg,jit_ci] = jitter_average(data,orig_times,max_jitter,maxlag)

num_reps = 1000;
t_jit_avg = zeros(num_reps,2*maxlag+1);
for r = 1:num_reps

    jit_calcs = zeros(length(orig_times),2*maxlag+1);
    for i = 1:length(orig_times)     
       
        if orig_times(i) > maxlag & orig_times(i) < length(data)-maxlag
            jit_val = ceil(rand*2*max_jitter)-max_jitter;
            jit_calcs(i,:) = data(orig_times(i)+jit_val-maxlag:orig_times(i)+jit_val+maxlag);
        
        else
            jit_calcs(i,:) = nan;
        end

    end

    t_jit_avg(r,:) = nanmean(jit_calcs);
    
end
jit_avg = mean(t_jit_avg);
t_jit_avg = sort(t_jit_avg);
jit_ci(1,:) = t_jit_avg(10,:);
jit_ci(2,:) = t_jit_avg(990,:);