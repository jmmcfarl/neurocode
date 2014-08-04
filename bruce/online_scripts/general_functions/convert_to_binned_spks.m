function binned_spikes = convert_to_binned_spks(spike_bins,NT)

T = tabulate(spike_bins);
binned_spikes = zeros(NT,1);
binned_spikes(T(:,1)) = T(:,2);