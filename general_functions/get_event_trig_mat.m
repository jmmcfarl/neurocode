function [ev_mat,lags] = get_event_trig_mat(sig,event_ids,backlag,forwardlag)

lags = -backlag:forwardlag;

[NT,p] = size(sig);

n_events = length(event_ids);
ev_mat = zeros(n_events,length(lags),p);

bad_ids = find(event_ids < backlag);
if ~isempty(bad_ids)
    fprintf('Dropping %d early events\n',length(bad_ids));
    event_ids(bad_ids) = [];
end
bad_ids = find(event_ids > length(sig) - forwardlag);
if ~isempty(bad_ids)
    fprintf('Dropping %d late events\n',length(bad_ids));
    event_ids(bad_ids) = [];
end
n_events = length(event_ids);

for i = 1:n_events
    %     fprintf('%d of %d\n',i,n_events);
    cur_ids = (event_ids(i)-backlag):(event_ids(i)+forwardlag);
    ev_mat(i,:,:) = sig(cur_ids,:);
end
