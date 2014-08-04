function store_spike_data(Spikes,dat_name,verbose)

if nargin < 3
    verbose = 0;
end

tic
Vrange = minmax(Spikes.V);
Spikes.V = int16((Spikes.V - Vrange(1))*2^16/diff(Vrange) - 2^15);
Spikes.Vrange = Vrange;
save(dat_name,'Spikes');
fprintf('Saving to %s, took %.3f sec\n',dat_name,toc);