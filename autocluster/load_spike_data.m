function Spikes = load_spike_data(dat_name)

load(dat_name);
Spikes.V = (double(Spikes.V) + 2^15)*diff(Spikes.Vrange)/2^16 + Spikes.Vrange(1);

