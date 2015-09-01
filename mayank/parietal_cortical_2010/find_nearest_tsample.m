function [t_samples] = find_nearest_tsample(t_axis,times)

t_samples = zeros(size(times));
for i = 1:length(times)
   
    [dummy,t_samples(i)] = min(abs(t_axis-times(i)));
    
end