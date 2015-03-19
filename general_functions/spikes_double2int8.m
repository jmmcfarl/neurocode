function spikes_out = spikes_double2int8(spikes_in)

nanvals = isnan(spikes_in);
spikes_out = int8(spikes_in);
spikes_out(nanvals) = 126;