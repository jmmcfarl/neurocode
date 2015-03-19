function spikes_out = spikes_int82double(spikes_in)

nanvals = find(spikes_in == 126);
spikes_out = double(spikes_in);
spikes_out(nanvals) = nan;
