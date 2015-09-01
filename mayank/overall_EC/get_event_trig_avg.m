function [event_trig_avg,data_mat] = get_event_trig_avg(data,event_samps,forwardlag,backwardlag)

if nargout > 1
    data_mat = nan(length(event_samps),forwardlag+backwardlag+1);
    for i = 1:length(event_samps)
        if event_samps(i) > backwardlag && event_samps(i) + forwardlag < length(data)
            data_mat(i,:) = data(event_samps(i)-backwardlag:event_samps(i)+forwardlag);
        end
    end
    event_trig_avg = nanmean(data_mat);
else
    event_trig_avg = 0;
    for i = 1:length(event_samps)
        if event_samps(i) > backwardlag && event_samps(i) + forwardlag < length(data)
            event_trig_avg = event_trig_avg + ...
                data(event_samps(i)-backwardlag:event_samps(i)+forwardlag);
        end
    end
    event_trig_avg = event_trig_avg/length(event_samps);
end