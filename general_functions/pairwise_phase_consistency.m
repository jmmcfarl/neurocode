function [ppc0,ppc1] = pairwise_phase_consistency(spike_phases,spike_trials)


max_spks = 10000;
n_spks = length(spike_phases);
if n_spks > max_spks
    error('No support for this many spikes')
end

spike_phases = spike_phases(:);

spike_phase_vecs = [cos(spike_phases) sin(spike_phases)];

spike_phase_dist = pdist(spike_phase_vecs,@phase_dist);
trial_diff = pdist(spike_trials(:),'hamming');
ppc0 = mean(spike_phase_dist);
ppc1 = mean(spike_phase_dist(trial_diff == 1));

end

function D = phase_dist(XI,XJ)
D = XJ(:,1)*XI(1) + XJ(:,2)*XI(2);
end
