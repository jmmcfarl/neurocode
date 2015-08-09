function [up_amps,down_amps] = get_state_amplitudes_mean(sig,up_inds,down_inds,trans_buffer)

if nargin < 3
    trans_buffer = 0;
end

up_amps = nan(size(up_inds));
down_amps = nan(size(down_inds));
n_states = length(up_inds);

for i = 1:n_states-1
    cur_inds = (up_inds(i) + trans_buffer):(down_inds(i)-trans_buffer); %set of indices in current state minus transition buffer
    up_amps(i) = mean(sig(cur_inds));
    cur_inds = (down_inds(i) + trans_buffer):(up_inds(i+1)-trans_buffer);
    down_amps(i) = mean(sig(cur_inds));
end
cur_inds = (up_inds(end) + trans_buffer):(down_inds(end)-trans_buffer); %set of indices in current state minus transition buffer
up_amps(end) = mean(sig(cur_inds));
