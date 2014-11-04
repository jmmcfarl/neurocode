function [up_amps,down_amps] = get_state_amplitudes(sig,up_inds,down_inds)

up_amps = nan(size(up_inds));
down_amps = nan(size(down_inds));
n_states = length(up_inds);

for i = 1:n_states-1
    up_amps(i) = median(sig(up_inds(i):down_inds(i)));
    down_amps(i) = median(sig(down_inds(i):up_inds(i+1)));
end
up_amps(end) = median(sig(up_inds(end):down_inds(end)));
